import m from "mithril";
import Chart from "chart.js/auto"
import md5 from "md5";

import Navbar from "../components/navbar";
import PtableSidebar from "../components/ptableSidebar";
import {
    updateData,
    fillMsDataGaps,
    smiDrawerTheme,
    getChartStyles,
    getHintColor,
    switchGraph,
    resizeContainer,
    initHint,
    initOptionHighlight,
    initMsChart,
    initSubmitBtn,
    initNextBtn,
    initMsQuestionSubmitBtn
} from "../libraries/home_helpers";
import {applyFilter, endTrial, initiateTutorial} from "../libraries/tutorial_helper";


export var Tutorial = {
    view: () => (
        <div className="content">
            <div id="overlay" class="black-overlay show">
                <div id="tutorial-text-box1" class="tutorial-text-box">
                    <h1 id="tutorial-text1" class="tutorial-text">Welcome to Hydragen</h1>
                </div>
                <div id="tutorial-text-box2" class="tutorial-text-box">
                    <h4 id="tutorial-text2" class="tutorial-text">Click Anywhere to Continue</h4>
                </div>
            </div>

            <div id="main-container">
                <Navbar/>
                <div className="container-wrapper">
                    <div className="container">
                        <form class="block is-relative" id="msQuestionForm">
                            <div className="column-left">
                                <div className="title-block">
                                    <h1>Mass Spectrometry Practice</h1>
                                </div>

                                <div className="chartPTableContainer">
                                    <div className="chartContainer">
                                        <p>Identify the compound that would give this mass spectrum.</p>
                                        <canvas id="msChart">
                                            <p>Loading...</p>
                                        </canvas>
                                    </div>
                                    <PtableSidebar/>
                                </div>

                                <button id="switch-graph" class="button is-primary">View Periodic Table</button>

                                <div className="feedback">
                                    <h4 id="question-feedback"></h4>
                                </div>
                            </div>
                            <div className="column-right">
                                <div className="answer-topbar">
                                    <div className="hint-container">
                                        <div id="hint-display" className="hint-display">
                                            <strong class="green-text" id="hints-danger-level"><span
                                                id="hints-used">3</span>/3 Hints
                                                Left</strong>
                                            <span class="tooltip" id="tooltip-text">Click on a bar in the MS chart to receive a hint about the fragments it usually represents.
                                        Note that not all peaks will have a corresponding hint, and that the hints are not necessarily correct for every molecule.</span>
                                        </div>

                                    </div>
                                    <div className="skip-container">
                                        <button class="move-on skip button" type="button" id="next"></button>
                                    </div>
                                </div>
                                <div id="answer" className="answer">
                                    <div className="control column mcq-options">
                                        <label class="radio">
                                            <input type="radio" name="answer" value="0"></input>
                                            <img id="radio-opt0" class="options"></img>
                                        </label>
                                    </div>
                                    <div className="control column mcq-options">
                                        <label class="radio">
                                            <input type="radio" name="answer" value="1"></input>
                                            <img id="radio-opt1" class="options"></img>
                                        </label>
                                    </div>
                                    <div className="control column mcq-options">
                                        <label class="radio">
                                            <input type="radio" name="answer" value="2"></input>
                                            <img id="radio-opt2" class="options"></img>
                                        </label>
                                    </div>
                                    <div className="control column mcq-options">
                                        <label class="radio">
                                            <input type="radio" name="answer" value="3"></input>
                                            <img id="radio-opt3" class="options"></img>
                                        </label>
                                    </div>
                                </div>
                                <div className="control submit-container">
                                    <input class="button is-primary" type="submit" id="submit"
                                           value="Submit Answer"></input>
                                </div>
                            </div>
                            <div className="is-overlay mcq-overlay is-hidden">
                                <h2>Loading...</h2>
                            </div>
                        </form>
                    </div>
                </div>
            </div>

        </div>
    ),
    oncreate: async () => {
        $("#switch-graph").on("click", async () => {
                switchGraph();
            }
        );
        //Respond to window size changes
        resizeContainer();
        window.addEventListener("resize", resizeContainer);

        initiateTutorial();

        // localStorage.setItem("tutorialQuestionsCompleted", 0); //reset no of trys (for testing only)
        if (localStorage.getItem("tutorialQuestionsCompleted") === null) {
            localStorage.setItem("tutorialQuestionsCompleted", 0);
        }
        ////End tutorial when user completed 3 questions (before)
        if (parseInt(localStorage.getItem("tutorialQuestionsCompleted")) >= 3) {
            console.log("end trial");
            endTrial();
        }

        await initHint();
        initOptionHighlight();
        await initMsChart();
        initMsQuestionSubmitBtn();
        initSubmitBtn(true);
        initNextBtn();
    }
}

export default Tutorial;
