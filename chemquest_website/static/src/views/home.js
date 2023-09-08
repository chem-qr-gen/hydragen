import m from "mithril";
import Chart from "chart.js/auto"
import md5 from "md5";

import Navbar from "../components/navbar";
import PtableSidebar from "../components/ptableSidebar";
import updateData from "../libraries/chartjs_helpers";

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

var MCQ = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <PtableSidebar />
                <div class="container">
                    <div class="block">
                        <h1>Mass Spectrometry Practice</h1>
                    </div>
                    <form class="block is-relative" id="msQuestionForm">
                        <input type="hidden" id="csrf_token"></input>
                        <div class="field">
                            <p>Identify the compound that would give this mass spectrum.</p>
                            <div class="columns">
                                <div class="column is-three-quarters">
                                    <canvas id="msChart">
                                        <p>Loading...</p>
                                    </canvas>
                                </div>
                                <div class="column is-one-quarter has-background-grey-lighter">
                                    <h5>Hints</h5>
                                    <p><i>Click on a bar in the MS chart to receive a hint about the fragments it usually represents.</i></p>
                                    <p><i>Note that not all peaks will have a corresponding hint, and that the hints are not necessarily correct for every molecule.</i></p>
                                    <p><i>You have used <b><span id="hints-used">0</span> of 3</b> hints.</i></p>
                                </div>
                            </div>
                        </div>
                        <div class="field">
                            <div class="columns">
                                <div class="control column">
                                    <label class="radio">
                                        <input type="radio" name="answer" value="0"></input>
                                        <img id="radio-opt0"></img>
                                    </label>
                                </div>
                                <div class="control column">
                                    <label class="radio">
                                        <input type="radio" name="answer" value="1"></input>
                                        <img id="radio-opt1"></img>
                                    </label>
                                </div>
                                <div class="control column">
                                    <label class="radio">
                                        <input type="radio" name="answer" value="2"></input>
                                        <img id="radio-opt2"></img>
                                    </label>
                                </div>
                                <div class="control column">
                                    <label class="radio">
                                        <input type="radio" name="answer" value="3"></input>
                                        <img id="radio-opt3"></img>
                                    </label>
                                </div>
                            </div>
                            <h4 id="question-feedback"></h4>
                        </div>
                        <div class="field is-grouped">
                            <div class="control">
                                <input class="button is-primary" type="submit" id="submit" value="Submit Answer"></input>
                            </div>
                            <div class="control">
                                <span class="button" id="next">Next Question</span>
                            </div>
                        </div>
                        <div class="is-overlay mcq-overlay is-hidden">
                            <h2>Loading...</h2>
                        </div>
                    </form>
                </div>
            </div>
        </div>
    ),
    oncreate: async () => {
        m.request({
            method: "GET",
            url: "/get_csrf_token"
        }).then(response => {
            $("#csrf_token").val(response.csrf_token);
        });

        // get hint data from server side
        var hintData = await m.request({
            method: "GET",
            url: "/get_hints"
        });
        var hintsUsed = 0;
        $("#hints-used").text(hintsUsed);

        // setup SMILES drawer
        var drawer = new SmiDrawer({
            themes: {
                dark: {
                    C: '#fff',
                    O: '#fff',
                    N: '#fff',
                    F: '#fff',
                    CL: '#fff',
                    BR: '#fff',
                    I: '#fff',
                    P: '#fff',
                    S: '#fff',
                    B: '#fff',
                    SI: '#fff',
                    H: '#fff',
                    BACKGROUND: '#141414'
                },
                light: {
                    C: '#222',
                    O: '#222',
                    N: '#222',
                    F: '#222',
                    CL: '#222',
                    BR: '#222',
                    I: '#222',
                    P: '#222',
                    S: '#222',
                    B: '#222',
                    SI: '#222',
                    H: '#222',
                    BACKGROUND: '#fff'
                }
            }
        });

        // highlight the selected option when clicked
        $("input[name='answer']").on("click", () => {
            $("input[name='answer']:checked").parent().parent()
            .addClass("radio-selected")
            .siblings().removeClass("radio-selected");
        });

        const getNewQuestion = async () => {
            // get new question with SMILES, mass spec data, etc.
            var data = await m.request({
                method: "GET",
                url: "/ms_questions_new",
                params: {id: "random"}
            });

            // filled data for the mass spec graph
            var filledMsData = fillMsDataGaps(data['ms_data']);

            // get random MCQ options
            var mcqAnswers = await m.request({
                method: "GET",
                url: "/generate_mcq",
                params: {input_smiles: data["smiles"]}
            }).then(response => response.map(i => i["SMILES"]));

            // insert the correct answer randomly into the options
            var correctAnswer = Math.floor(Math.random() * 3); 
            mcqAnswers.splice(correctAnswer, 0, data["smiles"]);

            // generate the structure images for display based on the SMILES MCQ options
            for (var index = 0; index < mcqAnswers.length; index++) {
                drawer.draw(mcqAnswers[index], "#radio-opt" + index);
            }

            // create a unique hash id from the timestamp and SMILES of the correct answer
            var timestamp = new Date().getTime();
            var questionHash = md5(timestamp + data["smiles"]);

            return {
                hashId: questionHash,
                rawData: data,
                filledMsData: filledMsData,
                mcqAnswers: mcqAnswers,
                correctAnswer: correctAnswer
            }
        };

        var questionData = await getNewQuestion();

        // chart initialisation
        var msChart = new Chart(
            document.getElementById("msChart"),
            {
                type: "bar",
                data: {
                    labels: questionData.filledMsData.map(entry => entry.mz),
                    datasets: [{
                        label: "Relative Abundance",
                        data: questionData.filledMsData.map(entry => entry.abundance),
                        backgroundColor: Array(questionData.filledMsData.length).fill("#000000"),
                        borderColor: Array(questionData.filledMsData.length).fill("#000000"),
                        barPercentage: 0.5
                    }]
                },
                options: {
                    plugins: {
                        legend: {
                            display: false
                        },
                        tooltip: {
                            titleFont: {
                                size: 16
                            },
                            bodyFont: {
                                size: 16
                            },
                            callbacks: {
                                label: (context) => {
                                    var mz = parseInt(context.label);
                                    var main_text = "Relative Abundance: " + context.parsed.y;
                                    if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == "#009900") {
                                        var hint = hintData.find(i => i.mz == mz);
                                        var hint_text = hint.hint_text;
                                    } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == "#990000") {
                                        var hint_text = "No hint available for this value."
                                    } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == "#996600") {
                                        var hint_text = "You have used all of your hints."
                                    } else {
                                        var hint_text = "Click to see a hint."
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
                                msChart.data.datasets[0].backgroundColor[index] = "#996600";
                            } else if (hint) {
                                msChart.data.datasets[0].backgroundColor[index] = "#009900";
                                hintsUsed++;
                                $("#hints-used").text(hintsUsed);
                                
                            } else {
                                msChart.data.datasets[0].backgroundColor[index] = "#990000";
                            }

                            msChart.update();
                        }
                    }
                }
            }
        );
        
        $("#next").on("click", async () => {
            // remove disabled from submit button
            $("#submit").prop("disabled", false);

            // reset the option buttons
            $("input[name='answer']").each((index, element) => {
                $(element).prop("checked", false);
                $(element).parent().parent().removeClass("radio-selected");
            });

            // reset the feedback text
            $("#question-feedback").removeClass("is-success is-danger");
            $("#question-feedback").text("");

            // reset hint data
            msChart.data.datasets[0].backgroundColor = ['rgba(0, 0, 0, 1)'];
            hintsUsed = 0;
            $("#hints-used").text(hintsUsed);

            // show loading overlay
            $(".mcq-overlay").removeClass("is-hidden");

            // get new question data and remove loading overlay
            questionData = await getNewQuestion();
            $(".mcq-overlay").addClass("is-hidden");
            updateData(msChart, questionData.filledMsData);
        });

        var isCorrect;

        $("#msQuestionForm").on("submit", () => {
            if (questionData.correctAnswer == $("input[name='answer']:checked").val()) {
                $("#question-feedback").removeClass("is-danger");
                $("#question-feedback").addClass("is-success");
                $("#question-feedback").text("Correct! Well done!");
                isCorrect = true;
            } else {
                $("#question-feedback").removeClass("is-success");
                $("#question-feedback").addClass("is-danger");
                $("#question-feedback").text("Incorrect, try again.");
                isCorrect = false;
            }

            // disable the submit button if the answer is correct
            if (isCorrect) {
                $("#submit").prop("disabled", true);
            }

            // submits the attempt to the server. Elo calculation and compiling of attempts for the same question are done server-side
            // TODO: add number of hints used
            m.request({
                method: "POST",
                url: "/record_attempt",
                headers: localStorage.getItem("jwt") ? {"Authorization": "Bearer " + localStorage.getItem("jwt")} : {},
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "hashId": questionData.hashId,
                    "qid": questionData.rawData.qid,
                    "options": questionData.mcqAnswers,
                    "answer": $("input[name='answer']:checked").val(),
                    "isCorrect": isCorrect
                }
            });
            return false;
        });
    }
}

export default MCQ;