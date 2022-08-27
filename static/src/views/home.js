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
        <div class="content container">
            <div class="block">
                <h1>Mass Spectrometry Practice</h1>
            </div>
            <div class="block">
                <div class="field">
                    <p id="question-id" style="display:none;"></p>
                    <p id="question-text-1">test 1</p>
                    <img id="question-img" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEAAAAALAAAAAABAAEAAAIBAAA="></img>
                    <p id="question-text-2">test 2</p>
                </div>
                <div class="field">
                    <div class="control">
                        <input class="input" type="text" placeholder="Answer" id="answer"></input>
                    </div>
                    <p id="question-feedback"></p>
                </div>
                <div class="field is-grouped">
                    <div class="control">
                        <button class="button is-primary" id="submit">Submit</button>
                    </div>
                    <div class="control">
                        <button class="button" id="next">Next Question</button>
                    </div>
                </div>
            </div>
            <div class="block tabs-with-content">
                <div class="tabs">
                    <ul>
                        <li><a>Hint 1</a></li>
                        <li><a>Hint 2</a></li>
                    </ul>
                </div>
                <section class="tab-content"></section>
                <section class="tab-content"></section>
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
                $("#answer").removeClass("is-success is-danger");
                $("#answer").val("");
                $("#question-feedback").removeClass("is-success is-danger");
                $("#question-id").text(response.id);
                $("#question-text-1").text(response.text1);
                $("#question-img").attr("src", response.imgsrc);
                $("#question-text-2").text(response.text2);
                answers = response.answers.split(";");
            });
        });
        $("#submit").on("click", () => {
            if (answers.includes($("#answer").val())) {
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

        var tabs = $(".tabs li")
        var tabsContent = $(".tab-content")

        const deactivateTabs = () => {
            tabs.each((index, tab) => {
                $(tab).removeClass("is-active");
            });
        }
        const hideTabsContent = () => {
            tabsContent.each((index, content) => {
                $(content).removeClass("is-active")
            });
        }
        tabs.each((index, tab) => {
            $(tab).on("click", () => {
                deactivateTabs();
                hideTabsContent();
                $(tab).addClass("is-active");
                tabsContent.eq(index).addClass("is-active");
            });
        });
    }
}

export default Home;
