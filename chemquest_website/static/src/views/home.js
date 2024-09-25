import m from "mithril";
import Chart from "chart.js/auto"
import md5 from "md5";
import JSConfetti from 'js-confetti';

import Navbar from "../components/navbar";
import PtableSidebar from "../components/ptableSidebar";
import {
  updateData,
  fillMsDataGaps,
  smiDrawerTheme,
  getChartStyles,
  getHintColor,
  switchGraph, resizeContainer
} from "../libraries/home_helpers";
import {applyFilter} from "../libraries/tutorial_helper";


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
                    <strong class="mschart-hint green-text" id="hints-danger-level"><span id="hints-used">3</span>/3 Hints
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
                  <input class="button is-primary" type="submit" id="submit" value="Submit Answer"></input>
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

    // get hint data from server side
    var hintData = await m.request({
      method: "GET",
      url: "/get_hints"
    });
    var hintsUsed = 0;
    $("#hints-used").text(3 - hintsUsed);
    getHintColor(hintsUsed);


    // setup SMILES drawer
    var drawer = new SmiDrawer({
      themes: smiDrawerTheme,
    });

    // highlight the selected option when clicked
    $("input[name='answer']").on("click", () => {
      $("input[name='answer']:checked").parent().parent()
        .addClass("radio-selected")
        .siblings().removeClass("radio-selected");
    });

    var generateTimestamp;
    const getNewQuestion = async () => {
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
        drawer.draw(mcqAnswers[index]["smiles"], "#radio-opt" + index);
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

    var questionData = await getNewQuestion();

    // chart initialisation
    const chartStyles = getChartStyles();
    var msChart = new Chart(
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
            },
            y: {
              grid: {color: chartStyles.gridColor},
              border: {color: chartStyles.axisColor},
              ticks: {color: chartStyles.tickColor}
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
                  if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.hintColor) {
                    var hint = hintData.find(i => i.mz == mz);
                    var hint_text = hint.hint_text;
                  } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.noHintColor) {
                    var hint_text = "No hint available for this value."
                  } else if (msChart.data.datasets[0].backgroundColor[context.dataIndex] == chartStyles.usedAllHintsColor) {
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
                msChart.data.datasets[0].backgroundColor[index] = chartStyles.usedAllHintsColor;
              } else if (hint) {
                if (msChart.data.datasets[0].backgroundColor[index] != chartStyles.hintColor){ //have not yet selected this element for hint
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

    const jsConfetti = new JSConfetti();
    var isCorrect;

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
          "isCorrect": isCorrect,
          "generate_timestamp": generateTimestamp
        }
      });
      return false;
    });

    $("#submit").on("click", async () => {
      if (questionData.correctAnswer == $("input[name='answer']:checked").val()) {    //Answer is correct
        $("#question-feedback").removeClass("is-danger");                           //Change Feedback
        $("#question-feedback").addClass("is-success");
        $("#question-feedback").text("Correct! Well done!");
        $(".move-on").addClass("next");                                             //Change skip button text
        $(".move-on").removeClass("skip")
        isCorrect = true;
        jsConfetti.addConfetti({
          confettiColors: ["#ea9a0c", "#2dc5f6", "#6046fe"]
        });
      } else {                                                                        //Answer Incorrect
        $("#question-feedback").removeClass("is-success");                          //Change Feedback
        $("#question-feedback").addClass("is-danger");
        $("#question-feedback").text("Incorrect, try again.");
        isCorrect = false;
        $("#msQuestionForm").shake();
      }

      // disable the submit button if the answer is correct
      if (isCorrect) {
        $("#submit").prop("disabled", true);
      }
    });

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
      console.log("removing text");
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
}

export default MCQ;
