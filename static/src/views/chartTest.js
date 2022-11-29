import m from "mithril";
import Chart from "chart.js/auto"

import Navbar from "../components/navbar";

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
                    <canvas id="msChart">
                        <p>Loading...</p>
                    </canvas>
                </div>
            </div>
        </div>
    ),
    oncreate: async () => {
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
    }
}

export default ChartTest;