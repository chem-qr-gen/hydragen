// home site template, reused in tutorial
import m from "mithril";
import Chart from "chart.js/auto";
import JSConfetti from "js-confetti";
import md5 from "md5";
import {endTrial} from "./tutorial_helper";
import {getHintColor, getChartStyles} from "./home_helpers"
import {getTheme} from "./universal_helpers";
import { smiDrawerTheme } from "./themes";

//Declare global variables
let hintsUsed;
let hintData
let questionData;
let isCorrect;
let generateTimestamp;
let msChart;
const chartStyles = getChartStyles();

async function initHint() {
    // get hint data from server side
    hintData = await m.request({
        method: "GET",
        url: "/get_hints"
    });
    hintsUsed = 0;
    $("#hints-used").text(3 - hintsUsed);
    getHintColor(hintsUsed);
}

function initOptionHighlight() {
    // highlight the selected option when clicked
    $("input[name='answer']").on("click", () => {
        $("input[name='answer']:checked").parent().parent()
            .addClass("radio-selected")
            .siblings().removeClass("radio-selected");
    });
}

const getNewQuestion = async () => {
    var theme = getTheme();
    // setup SMILES drawer
    var drawer = new SmiDrawer({
        themes: smiDrawerTheme
    });
    // get new question with SMILES, mass spec data, etc.
    var data = await m.request({
        method: "GET",
        url: "/ms_questions_new",
        params: {id: "random"}
    });

    // get generation time (to be passed for record_attempt)
    generateTimestamp = data['generate_timestamp'];

    // filled data for the mass spec graph
    var filledMsData = fillMsDataGaps(data['ms_data']);

    // get random MCQ options
    var mcqAnswers = await m.request({
        method: "GET",
        url: "/generate_mcq",
        params: {input_smiles: data["smiles"]}
    });

    // insert the correct answer randomly into the options
    var correctAnswerIndex = Math.floor(Math.random() * 3);
    mcqAnswers.splice(correctAnswerIndex, 0, {
        "smiles": data["smiles"],
        "explanation": "correct"
    });

    // generate the structure images for display based on the SMILES MCQ options
    for (var index = 0; index < mcqAnswers.length; index++) {
        drawer.draw(mcqAnswers[index]["smiles"], "#radio-opt" + index, theme);
    }

    // create a unique hash id from the timestamp and SMILES of the correct answer
    var timestamp = new Date().getTime();
    var questionHash = md5(timestamp + data["smiles"]);

    return {
        hashId: questionHash,
        rawData: data,
        filledMsData: filledMsData,
        mcqAnswers: mcqAnswers,
        correctAnswer: correctAnswerIndex
    }
};

async function initMsChart() {
    questionData = await getNewQuestion();
    // chart initialisation
    msChart = new Chart(
        document.getElementById("msChart"),
        {
            type: "bar",
            data: {
                labels: questionData.filledMsData.map(entry => entry.mz),
                datasets: [{
                    label: "Relative Abundance",
                    data: questionData.filledMsData.map(entry => entry.abundance),
                    backgroundColor: Array(questionData.filledMsData.length).fill(chartStyles.backgroundColor),
                    borderColor: Array(questionData.filledMsData.length).fill(chartStyles.axisColor),
                    barPercentage: chartStyles.barPercentage
                }]
            },
            options: {
                scales: {
                    x: {
                        grid: {color: chartStyles.gridColor},
                        border: {color: chartStyles.axisColor},
                        ticks: {color: chartStyles.tickColor},
                        title: {
                            display: true,
                            text: "m/z",
                            color: chartStyles.labelColor,
                            font: {
                                family: chartStyles.labelFontFamily,
                                size: chartStyles.labelFontSize,
                                weight: chartStyles.labelFontWeight,
                            },
                        },
                    },
                    y: {
                        grid: {color: chartStyles.gridColor},
                        border: {color: chartStyles.axisColor},
                        ticks: {color: chartStyles.tickColor},
                        title: {
                            display: true,
                            text: "Relative Intensity",
                            color: chartStyles.labelColor,
                            font: {
                                family: chartStyles.labelFontFamily,
                                size: chartStyles.labelFontSize,
                                weight: chartStyles.labelFontWeight,
                            },
                        }
                    }
                },
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        titleFont: {
                            size: chartStyles.titleFontSize
                        },
                        titleColor: chartStyles.titleColor,
                        bodyFont: {
                            size: chartStyles.bodyFontSize
                        },
                        bodyColor: chartStyles.bodyColor,
                        callbacks: {
                            label: (context) => {
                                var mz = parseInt(context.label);
                                var main_text = "Relative Abundance: " + context.parsed.y;
                                var hint_text;
                                if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.hintColor) {
                                    var hint = hintData.find(i => i.mz == mz);
                                    hint_text = hint.hint_text;
                                } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.noHintColor) {
                                    hint_text = "No hint available for this value."
                                } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.usedAllHintsColor) {
                                    hint_text = "You have used all of your hints."
                                } else {
                                    hint_text = "Click to see a hint."
                                }
                                return [main_text, hint_text];
                            }
                        }
                    }
                },
                onClick: (event, elements) => {
                    if (elements.length > 0) {
                        var index = elements[0].index;

                        var mz = parseInt(msChart.data.labels[index]);
                        var hint = hintData.find(i => i.mz == mz);

                        if (hintsUsed >= 3) {
                            msChart.data.datasets[0].backgroundColor[index] = chartStyles.usedAllHintsColor;
                        } else if (hint) {
                            if (msChart.data.datasets[0].backgroundColor[index] != chartStyles.hintColor) { //have not yet selected this element for hint
                                msChart.data.datasets[0].backgroundColor[index] = chartStyles.hintColor;
                                hintsUsed++;
                                $("#hints-used").text(3 - hintsUsed);
                                getHintColor(hintsUsed);
                            }
                        } else {
                            msChart.data.datasets[0].backgroundColor[index] = chartStyles.noHintColor;
                        }

                        msChart.update();
                    }
                }
            }
        }
    );
}

function initMsQuestionSubmitBtn() {
    $("#msQuestionForm").on("submit", () => {
        // submits the attempt to the server. Elo calculation and compiling of attempts for the same question are done server-side
        // TODO: add number of hints used
        m.request({
            method: "POST",
            url: "/record_attempt",
            body: {
                "hashId": questionData.hashId,
                "qid": questionData.rawData.qid,
                "options": questionData.mcqAnswers,
                "answer": $("input[name='answer']:checked").val(),
                "hintsUsed": hintsUsed,
                "isCorrect": isCorrect,
                "generate_timestamp": generateTimestamp
            }
        });
        if (isCorrect) {
            $("#submit").prop("disabled", true); // disable the submit button
        }
        return false;
    });
}

function initSubmitBtn(isTutorial) {
    const jsConfetti = new JSConfetti();
    $("#submit").on("click", async () => {
        if (questionData.correctAnswer == $("input[name='answer']:checked").val()) { //Answer is correct
            $("#question-feedback").removeClass("is-danger"); //Change Feedback
            $("#question-feedback").addClass("is-success");
            $("#question-feedback").text("Correct! Well done!");
            $(".move-on").addClass("next"); //Change skip button text
            $(".move-on").removeClass("skip")
            isCorrect = true;
            await jsConfetti.addConfetti({
                confettiColors: ["#ea9a0c", "#2dc5f6", "#6046fe"]
            });
            if (isTutorial) {
                localStorage.setItem("tutorialQuestionsCompleted", parseInt(localStorage.getItem("tutorialQuestionsCompleted")) + 1) //Add 1 completed questions
                if (parseInt(localStorage.getItem("tutorialQuestionsCompleted")) >= 3) { //End tutorial when user completed 3 questions
                    endTrial();
                }
            }
        } else {//Answer Incorrect
            $("#question-feedback").removeClass("is-success");//Change Feedback
            $("#question-feedback").addClass("is-danger");

            // get the explanation tag for the chosen option
            var explan = questionData.mcqAnswers[$("input[name='answer']:checked").val()]["explanation"];

            // feedback based on type of error
            if (explan.some(i => ["halide-alcohol", "halide-amine", "halide-thiol", "alcohol-halide", "amine-halide", "thiol-halide"].includes(i))) { // halides - look at M+2/M+4 peaks 
            $("#question-feedback").text("Incorrect. Examine the peaks around the molecular ion.");
            } else if (explan.some(i => ["primary-o-flip", "secondary-amine-flip", "carbonyl-methyl", "ether-alcohol", "ester-flip", "ester-ketone-alcohol"].includes(i))) { // MW is the same, but the structure is different
            $("#question-feedback").text("Incorrect. Pay attention to the structure and fragment peaks.")
            } else if (explan.some(i => ["alcohol-amine", "alcohol-ether", "amine-alcohol", "secondary-amine-ether", "ester-amide", "acid-ester", "acid-amide", "nitrile-alkyne"].includes(i))) { // MW is different
            $("#question-feedback").text("Incorrect. Pay attention to the molecular ion peak.")
            } else if (explan.includes("similarity-algo")) { // similarity algorithm
            $("#question-feedback").text("Incorrect, but the correct structure is similar.")
            } else { // fallback
            $("#question-feedback").text("Incorrect, try again.");
            }
            isCorrect = false;
            $("#msQuestionForm").shake();
        }
    });
}

function initNextBtn() {
    $("#next").on("click", async () => {
        // remove disabled from submit button
        $("#submit").prop("disabled", false);

        // reset the option buttons
        $("input[name='answer']").each((index, element) => {
            $(element).prop("checked", false);
            $(element).parent().parent().removeClass("radio-selected");
        });
        // reset skip button
        $("#next").removeClass("next");
        $("#next").addClass("skip");

        // reset the feedback text
        $("#question-feedback").removeClass("is-success is-danger");
        $("#question-feedback").text("a");

        // reset hint data
        msChart.data.datasets[0].backgroundColor = [chartStyles.backgroundColor];
        hintsUsed = 0;
        $("#hints-used").text(3 - hintsUsed);
        getHintColor(hintsUsed);

        // show loading overlay
        $(".mcq-overlay").removeClass("is-hidden");

        // get new question data and remove loading overlay
        questionData = await getNewQuestion();
        $(".mcq-overlay").addClass("is-hidden");
        updateData(msChart, questionData.filledMsData);
    });
}

// update the data in a chart
const updateData = (chart, data) => {
    chart.data.labels = data.map(entry => entry.mz);
    chart.data.datasets[0].data = data.map(entry => entry.abundance);
    chart.update();
}

// fills "gaps" in the chart with zeros so it looks more like a proper MS chart
const fillMsDataGaps = msData => {
    // get the highest mz value in the data
    var highestMz = msData[msData.length - 1][0]
    var newMsData = []

    // fill the gaps with zeros
    for (var i = 1; i <= highestMz; i++) {
        newMsData.push({"mz": i, "abundance": 0});
    }
    // fill in the data from the server
    for (const i of msData) {
        newMsData[i[0] - 1]["abundance"] = i[1];
    }
    return newMsData;
}

export {
    initHint,
    initOptionHighlight,
    initMsChart,
    initMsQuestionSubmitBtn,
    initSubmitBtn,
    initNextBtn,
    updateData,
    fillMsDataGaps,
};