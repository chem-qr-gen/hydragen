import {getChartStylesFromTheme} from "./themes";
import {getTheme} from "./universal_helpers";

function getChartStyles() {
    var chartStyles = new Object();
    var theme = getTheme();
    chartStyles = getChartStylesFromTheme(chartStyles, theme);
    return chartStyles;
}

function switchGraph() {
    $(".ptable").toggleClass("hide");
    $(".chartContainer").toggleClass("hide");

    // change text of the button to "Graph" if periodic table is active, and "Periodic Table" if graph is active
    if ($(".ptable").hasClass("hide")) {
        $("#switch-graph").text("View Periodic Table");
    } else {
        $("#switch-graph").text("Back to Graph");
    }
}

function resizeContainer() {
    $("#main-container").css("width", document.documentElement.clientWidth)
    $("#main-container").css("height", document.documentElement.clientHeight)
}

// Update Hint indicator color
function getHintColor(hintsUsed) {
    switch (hintsUsed) {
        case (3 || 2): {
            $("#hints-danger-level").removeClass("yellow-text");
            $("#hints-danger-level").removeClass("green-text");
            $("#hints-danger-level").addClass("red-text");
        }
            break;
        case 1: {
            $("#hints-danger-level").removeClass("red-text");
            $("#hints-danger-level").removeClass("green-text");
            $("#hints-danger-level").addClass("yellow-text")
        }
            break;
        case 0: {
            $("#hints-danger-level").removeClass("yellow-text");
            $("#hints-danger-level").removeClass("red-text");
            $("#hints-danger-level").addClass("green-text");
        }
    }
    return;
}

// jQuery shake effect for wrong answer
// source: https://stackoverflow.com/questions/4399005/implementing-jquerys-shake-effect-with-animate
$.fn.shake = function (interval, distance, times) {
    interval = typeof interval == "undefined" ? 100 : interval;
    distance = typeof distance == "undefined" ? 10 : distance;
    times = typeof times == "undefined" ? 3 : times;
    var jTarget = $(this);
    jTarget.css('position', 'relative');
    for (var iter = 0; iter < (times + 1); iter++) {
        jTarget.animate({left: ((iter % 2 == 0 ? distance : distance * -1))}, interval);
    }
    return jTarget.animate({left: 0}, interval);
}

export {
    switchGraph,
    resizeContainer,
    getChartStyles,
    getHintColor,
};

