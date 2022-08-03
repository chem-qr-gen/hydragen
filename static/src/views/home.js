import m from "mithril";
import $ from "cash-dom";

var answers;
const getNewQuestion = () => {
    return m.request({
        method: "GET",
        url: "/ms_questions",
        params: {id: "random"}
    });
}

var Home = {
    view: () => (
        <div class="container">
            <h3>Mass Spectrometry Practice</h3>
            <div>
                <fieldset>
                    <p id="question-id" style="display:none;"></p>
                    <p id="question-text-1">test 1</p>
                    <img id="question-img" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEAAAAALAAAAAABAAEAAAIBAAA="></img>
                    <p id="question-text-2">test 2</p>
                    <input type="text" placeholder="Answer" id="answer"></input>
                    <p id="question-feedback"></p>
                    <button class="button" id="submit">Submit</button>
                    &nbsp;
                    <button class="button button-outline" id="next">Next Question</button>
                </fieldset>
            </div>
        </div>
    ),
    oncreate: () => {
        getNewQuestion().then(response => {
            $("#question-id").text(response.id);
            $("#question-text-1").text(response.text1);
            $("#question-img").attr("src", response.imgsrc);
            $("#question-text-2").text(response.text2);
            answers = response.answers.split(";");
        });
        $("#next").on("click", () => {
            getNewQuestion().then(response => {
                $("#answer").removeClass("wrong correct");
                $("#answer").val("");
                $("#question-feedback").removeClass("wrong correct");
                $("#question-id").text(response.id);
                $("#question-text-1").text(response.text1);
                $("#question-img").attr("src", response.imgsrc);
                $("#question-text-2").text(response.text2);
                answers = response.answers.split(";");
            });
        });
        $("#submit").on("click", () => {
            if (answers.includes($("#answer").val())) {
                $("#answer").removeClass("wrong");
                $("#answer").addClass("correct");
                $("#question-feedback").removeClass("wrong");
                $("#question-feedback").addClass("correct");
                $("#question-feedback").text("Correct! Well done!");
            } else {
                $("#answer").removeClass("correct");
                $("#answer").addClass("wrong");
                $("#question-feedback").removeClass("correct");
                $("#question-feedback").addClass("wrong");
                $("#question-feedback").text("Incorrect, try again.");
            }
        });
    }
}

export default Home;