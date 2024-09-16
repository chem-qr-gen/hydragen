//Define Chart Styles
//chartStyles: [
//    0 backgroundColor, 1 borderColor, 2 barPercentage, 3 titleFont size,
//    4 bodyFont size, 5 hint color, 6 noHint color, 7 usedAllHints color,
//    8 titleColor, 9 bodyColor, 10 gridColor, 11 axis color, 12 tick color
//]
function getChartStyles() {
  const chartStyles = new Object();
  const theme = document.documentElement.getAttribute('data-theme');
  console.log(theme);
  let _ = chartStyles;
  switch (theme) {
    case 'auto':
      if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
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
        _.gridColor = '#383B4A';
        _.axisColor = '#2a2a2a';
        _.tickColor = '#454545';

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
      _.gridColor = '#383B4A';
      _.axisColor = '#2a2a2a';
      _.tickColor = '#454545';
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

export {
  switchGraph,
  updateData,
  fillMsDataGaps,
  smiDrawerTheme,
  getChartStyles,
  getHintColor,
};

