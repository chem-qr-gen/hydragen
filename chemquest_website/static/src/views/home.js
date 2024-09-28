import m from "mithril";
import Chart from "chart.js/auto"
import md5 from "md5";
import JSConfetti from 'js-confetti';

import Navbar from "../components/navbar";
import PtableSidebar from "../components/ptableSidebar";
import {
    initNextBtn,
    initSubmitBtn,
    initMsChart,
    initHint,
    initOptionHighlight,
    initMsQuestionSubmitBtn
} from "../libraries/home_template";
import {switchGraph, resizeContainer } from "../libraries/home_helpers";

export var MCQ = {
    view: () => (
        <div className="content">
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

                                <div className="hint-container" style="width:100%; padding-bottom: 8px;">
                                    <div className="hint-display">
                                        <strong class="mschart-hint green-text" id="hints-danger-level"><span
                                            id="hints-used">3</span>/3
                                            Hints
                                            Left</strong>
                                        <span class="tooltip" id="tooltip-text">Click on a bar in the MS chart to receive a hint about the fragments it usually represents.
                                      Note that not all peaks will have a corresponding hint, and that the hints are not necessarily correct for every molecule.</span>
                                    </div>
                                </div>

                                <button id="switch-graph" class="button is-primary">View Periodic Table</button>

                                <div className="feedback">
                                    <h4 id="question-feedback"></h4>
                                </div>
                            </div>
                            <div className="column-right">
                                <div className="answer-topbar">

                                    <div className="skip-container">
                                        <button class="move-on skip button" type="button" id="next"></button>
                                    </div>
                                </div>
                                <div className="answer">
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

        await initHint();
        initOptionHighlight();
        await initMsChart();
        initMsQuestionSubmitBtn();
        initSubmitBtn(false);
        initNextBtn();
        
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

    }
}

export default MCQ;
