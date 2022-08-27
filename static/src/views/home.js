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
        var tabs = $(".tabs li")
        var tabsContent = $(".tab-content")
        var id = $("#question-id");
        var text1 = $("#question-text-1");
        var img = $("#question-img");
        var text2 = $("#question-text-2");
        var ansInput = $("#answer");
        var feedback = $("#question-feedback");

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

        const displayNewQuestion = () => {
            return getNewQuestion().then(response => {
                id.text(response.id);
                text1.text(response.text1);
                img.attr("src", response.imgsrc);
                text2.text(response.text2);
                tabsContent.eq(0).text(response.hint1);
                tabsContent.eq(1).text(response.hint2);
                answers = response.answers.split(";");
            });
        }

        displayNewQuestion();

        $("#next").on("click", () => {
            displayNewQuestion().then(() => {
                deactivateTabs();
                hideTabsContent();
                ansInput.removeClass("is-success is-danger");
                ansInput.val("");
                feedback.removeClass("is-success is-danger");
            });
        });

        $("#submit").on("click", () => {
            if (answers.includes(ansInput.val())) {
                ansInput.removeClass("is-danger");
                ansInput.addClass("is-success");
                feedback.removeClass("is-danger");
                feedback.addClass("is-success");
                feedback.text("Correct! Well done!");
            } else {
                ansInput.removeClass("is-success");
                ansInput.addClass("is-danger");
                feedback.removeClass("is-success");
                feedback.addClass("is-danger");
                feedback.text("Incorrect, try again.");
            }
        });
    }
}

export default Home;
