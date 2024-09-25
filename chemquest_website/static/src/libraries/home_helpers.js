//Define Chart Styles

function getChartStyles() {
  const chartStyles = new Object();
  const theme = document.documentElement.getAttribute('data-theme');
  console.log(theme);
  let _ = chartStyles;
  switch (theme) {
    case 'auto':
      if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
        _.backgroundColor = '#1a1752';
        _.borderColor = '#3f3f3f';
        _.barPeercentage = 0.5;
        _.titleFontSize = 16;
        _.titleColor = '#9da5cc';
        _.bodyFontSize = 16;
        _.bodyColor = '#9da5cc';
        _.hintColor = '#009900';
        _.noHintColor = '#990000';
        _.usedAllHintsColor = '#996600';
        _.gridColor = '#4247618f';
        _.axisColor = '#2a2a2a';
        _.tickColor = '#454545';
        _.labelColor = '#312000ff'
        _.labelFontFamily = 'Montserrat'
        _.labelFontSize = '16'
        _.labelFontWeight = '600'

      } else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
        _.backgroundColor = '#b6b6b6';
        _.borderColor = '#b6b6b6';
        _.barPeercentage = 0.5;
        _.titleFontSize = 16;
        _.titleColor = '#9da5cc';
        _.bodyFontSize = 16;
        _.bodyColor = '#9da5cc';
        _.hintColor = '#009900';
        _.noHintColor = '#990000';
        _.usedAllHintsColor = '#996600';
        _.gridColor = '#454545';
        _.axisColor = '#666979';
        _.tickColor = '#c8c8c8';
        _.labelColor = '#eec77eff'
        _.labelFontFamily = 'Montserrat'
        _.labelFontSize = '16'
        _.labelFontWeight = '600'
      }
      break;
    case 'light':
      _.backgroundColor = '#3f3f3f';
      _.borderColor = '#3f3f3f';
      _.barPeercentage = 0.5;
      _.titleFontSize = 16;
      _.titleColor = '#9da5cc';
      _.bodyFontSize = 16;
      _.bodyColor = '#9da5cc';
      _.hintColor = '#009900';
      _.noHintColor = '#990000';
      _.usedAllHintsColor = '#996600';
      _.gridColor = '#4247618f';
      _.axisColor = '#2a2a2a';
      _.tickColor = '#454545';
      _.labelColor = '#312000ff'
      _.labelFontFamily = 'Montserrat'
      _.labelFontSize = '16'
      _.labelFontWeight = '600'
      break;
    case 'dark':
      _.backgroundColor = '#b6b6b6';
      _.borderColor = '#b6b6b6';
      _.barPeercentage = 0.5;
      _.titleFontSize = 16;
      _.titleColor = '#9da5cc';
      _.bodyFontSize = 16;
      _.bodyColor = '#9da5cc';
      _.hintColor = '#009900';
      _.noHintColor = '#990000';
      _.usedAllHintsColor = '#996600';
      _.gridColor = '#454545';
      _.axisColor = '#666979';
      _.tickColor = '#c8c8c8';
      _.labelColor = '#eec77eff'
      _.labelFontFamily = 'Montserrat'
      _.labelFontSize = '16'
      _.labelFontWeight = '600'
      break;
  }
  console.log(chartStyles);
  return chartStyles
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
  console.log(hintsUsed);
  console.log($("#hints-danger-level"));
  switch (hintsUsed) {
    case (3 || 2): {
      $("#hints-danger-level").removeClass("yellow-text");
      $("#hints-danger-level").removeClass("green-text");
      $("#hints-danger-level").addClass("red-text");
    }
      console.log(2);
      break;
    case 1: {
      $("#hints-danger-level").removeClass("red-text");
      $("#hints-danger-level").removeClass("green-text");
      $("#hints-danger-level").addClass("yellow-text")
    }
      console.log(1);
      break;
    case 0: {
      $("#hints-danger-level").removeClass("yellow-text");
      $("#hints-danger-level").removeClass("red-text");
      $("#hints-danger-level").addClass("green-text");
    }
      console.log(0);
  }
  return;
}

// update the data in a chart
const updateData = (chart, data) => {
  chart.data.labels = data.map(entry => entry.mz);
  chart.data.datasets[0].data = data.map(entry => entry.abundance);
  chart.update();
}

// fills "gaps" in the chart with zeros so it looks more like a proper MS chart
const fillMsDataGaps = msData => {
  // get the highest mz value in the data
  var highestMz = msData[msData.length - 1][0]
  var newMsData = []

  // fill the gaps with zeros
  for (var i = 1; i <= highestMz; i++) {
    newMsData.push({"mz": i, "abundance": 0});
  }
  // fill in the data from the server
  for (const i of msData) {
    newMsData[i[0] - 1]["abundance"] = i[1];
  }
  return newMsData;
}

var smiDrawerTheme = {
  dark: {
    C: '#fff',
    O: '#fff',
    N: '#fff',
    F: '#fff',
    CL: '#fff',
    BR: '#fff',
    I: '#fff',
    P: '#fff',
    S: '#fff',
    B: '#fff',
    SI: '#fff',
    H: '#fff',
    BACKGROUND: '#141414'
  },
  light: {
    C: '#222',
    O: '#222',
    N: '#222',
    F: '#222',
    CL: '#222',
    BR: '#222',
    I: '#222',
    P: '#222',
    S: '#222',
    B: '#222',
    SI: '#222',
    H: '#222',
    BACKGROUND: '#fff'
  }
}

// jQuery shake effect for wrong answer
// source: https://stackoverflow.com/questions/4399005/implementing-jquerys-shake-effect-with-animate
$.fn.shake = function(interval,distance,times) {
  interval = typeof interval == "undefined" ? 100 : interval;
  distance = typeof distance == "undefined" ? 10 : distance;
  times = typeof times == "undefined" ? 3 : times;
  var jTarget = $(this);
  jTarget.css('position','relative');
  for(var iter=0;iter<(times+1);iter++){
     jTarget.animate({ left: ((iter%2==0 ? distance : distance*-1))}, interval);
  }
  return jTarget.animate({ left: 0},interval);
}

export {
  switchGraph,
  resizeContainer,
  updateData,
  fillMsDataGaps,
  smiDrawerTheme,
  getChartStyles,
  getHintColor,
};

