import m from "mithril";
import Chart from "chart.js/auto"

import Navbar from "../components/navbar";
import updateData from "../libraries/chartjs_helpers";

// fills "gaps" in the chart with zeros so it looks more like a proper MS chart
const fillMsDataGaps = msData => {
    var dataLength = msData[msData.length - 1]['mz']
    var newMsData = []
    for (var i = 1; i <= dataLength; i++) {
        newMsData.push({"mz": i, "abundance": 0});
    }
    for (const i of msData) {
        newMsData[i["mz"] - 1]["abundance"] = i["abundance"];
    }
    return newMsData;
}

var MCQ = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container">
                <div class="block">
                    <h1>Mass Spectrometry Practice</h1>
                </div>
                <form class="block" id="msQuestionForm">
                    <input type="hidden" id="csrf_token"></input>
                    <div class="field">
                        <p>Identify the compound that would give this mass spectrum.</p>
                        <div class="chartContainer">
                            <canvas id="msChart">
                                <p>Loading...</p>
                            </canvas>
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
                        <p id="question-feedback"></p>
                    </div>
                    <div class="field is-grouped">
                        <div class="control">
                            <input class="button is-primary" type="submit" id="submit" value="Submit Answer"></input>
                        </div>
                        <div class="control">
                            <span class="button" id="next">Next Question</span>
                        </div>
                    </div>
                </form>
                <div class="block">
                    
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

        var drawer = new SmiDrawer();

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
            var filledMsData = fillMsDataGaps(data['ms_data']);

            // get random MCQ options
            var mcqAnswers = await m.request({
                method: "GET",
                url: "/generate_mcq",
                params: {input_smiles: data["smiles"]}
            }).then(response => response.map(i => i["SMILES"]));

            var correctAnswer = Math.floor(Math.random() * 3); 
            mcqAnswers.splice(correctAnswer, 0, data["smiles"]); // insert the correct answer randomly into the options

            for (var index = 0; index < mcqAnswers.length; index++) {
                drawer.draw(mcqAnswers[index], "#radio-opt" + index); // generate the structure images based on the SMILES MCQ options
            }

            return {
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
                        backgroundColor: [
                            'rgba(0, 0, 0, 1)'
                        ],
                        borderColor: [
                            'rgba(0, 0, 0, 1)'
                        ],
                        barPercentage: 0.5
                    }]
                },
                options: {
                    plugins: {
                        legend: {
                            display: false
                        }
                    }
                }
            }
        );
        
        $("#next").on("click", async () => {
            $("input[name='answer']").each((index, element) => {
                $(element).prop("checked", false);
                $(element).parent().parent().removeClass("radio-selected");
            });
            $("#question-feedback").removeClass("is-success is-danger");
            $("#question-feedback").text("");

            questionData = await getNewQuestion();
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
            m.request({
                method: "POST",
                url: "/record_attempt",
                headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")},
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "id": questionData.rawData.qid,
                    "answer": questionData.mcqAnswers[$("input[name='answer']:checked").val()],
                    "isCorrect": isCorrect
                }
            });
            return false;
        });
    }
}

export default MCQ;