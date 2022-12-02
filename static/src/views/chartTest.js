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

var ChartTest = {
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
                        <canvas id="msChart">
                            <p>Loading...</p>
                        </canvas>
                    </div>
                    <div class="field">
                        <div class="control">
                            <input class="input" type="text" placeholder="Answer" id="answer" autocomplete="off"></input>
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

        var data = await m.request({
            method: "GET",
            url: "/ms_questions_new",
            params: {id: "random"}
        });
        var filledData = fillMsDataGaps(data['ms_data']);

        var msChart = new Chart(
            document.getElementById("msChart"),
            {
                type: "bar",
                data: {
                    labels: filledData.map(entry => entry.mz),
                    datasets: [{
                        label: "Relative Abundance",
                        data: filledData.map(entry => entry.abundance),
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
            data = await m.request({
                method: "GET",
                url: "/ms_questions_new",
                params: {id: "random"}
            });
            console.log(data.smiles);
            filledData = fillMsDataGaps(data['ms_data']);
            updateData(msChart, filledData);
        });

        document.getElementById("msQuestionForm").addEventListener("submit", e => {
            e.preventDefault();
            if ($("#answer").val() === data.smiles) {
                $("#answer").removeClass("is-danger");
                $("#answer").addClass("is-success");
                $("#question-feedback").removeClass("is-danger");
                $("#question-feedback").addClass("is-success");
                $("#question-feedback").text("Correct! Well done!");
            } else {
                $("#answer").removeClass("is-success");
                $("#answer").addClass("is-danger");
                $("#question-feedback").removeClass("is-success");
                $("#question-feedback").addClass("is-danger");
                $("#question-feedback").text("Incorrect, try again.");
            }
        });

    }
}

export default ChartTest;