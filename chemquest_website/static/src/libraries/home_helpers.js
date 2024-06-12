// Get Chart colors
//chartStyles structure: [
//    1 backgroundColor, 2 borderColor, 3 barPercentage, 4 titleFont size, 
//    5 bodyFont size, 6 hint color, 7 noHint color, 8 usedAllHints color,
//    9 titleColor, 10 bodyColor, 11 gridColor
//]
function getChartStyles() {
    const chartStyles = new Array(); 
    const theme = document.documentElement.getAttribute('data-theme');
    console.log(theme);
    switch(theme) {
        case 'auto':
            if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
                chartStyles.push('#000000','#000000', 0.5, 16, 16, '#009900', '#990000', "#996600", '#000000', '#000000', '#575d5e')
            }
            else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
                chartStyles.push('#FFFFFF','#FFFFFF', 0.5, 16, 16, '#009900', '#990000', "#996600",'#FFFFFF', '#FFFFFF', '#dbdbdb')
            }          
        case 'light':
            chartStyles.push('#000000','#000000', 0.5, 16, 16, '#009900', '#990000', "#996600", '#000000', '#000000', '#575d5e')
        case 'dark':
            chartStyles.push('#FFFFFF','#FFFFFF', 0.5, 16, 16, '#009900', '#990000', "#996600",'#FFFFFF', '#FFFFFF', '#dbdbdb')

    }
    console.log(chartStyles);
    return chartStyles
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
    updateData,
    fillMsDataGaps,
    smiDrawerTheme,
    getChartStyles
};

