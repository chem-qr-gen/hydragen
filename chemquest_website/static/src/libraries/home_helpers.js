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
    switch(theme) {
        case 'auto':
            if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
                _.backgroundColor = '#000000';
                _.borderColor = '#000000';
                _.barPeercentage = 0.5;
                _.titleFontSize = 16;
                _.titleColor = '#9da5cc';
                _.bodyFontSize = 16;
                _.bodyColor = '#9da5cc';
                _.hintColor = '#009900';
                _.noHintColor = '#990000';
                _.usedAllHintsColor = '#996600';
                _.gridColor = '#9da5cc';
                _.axisColor = '#212121';
                _.tickColor = '#212121';
            
            }
            else if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
                _.backgroundColor = '#FFFFFF';
                _.borderColor = '#FFFFFF';
                _.barPeercentage = 0.5;
                _.titleFontSize = 16;
                _.titleColor = '#9da5cc';
                _.bodyFontSize = 16;
                _.bodyColor = '#9da5cc';
                _.hintColor = '#009900';
                _.noHintColor = '#990000';
                _.usedAllHintsColor = '#996600';
                _.gridColor = '#9da5cc';
                _.axisColor = '#212121';
                _.tickColor = '#212121';        
            }  
        case 'light':
            _.backgroundColor = '#000000';
            _.borderColor = '#000000';
            _.barPeercentage = 0.5;
            _.titleFontSize = 16;
            _.titleColor = '#9da5cc';
            _.bodyFontSize = 16;
            _.bodyColor = '#9da5cc';
            _.hintColor = '#009900';
            _.noHintColor = '#990000';
            _.usedAllHintsColor = '#996600';
            _.gridColor = '#383B4A';
            _.axisColor = '#F5F5F5';
            _.tickColor = '#F5F5F5';
        
        case 'dark':
            _.backgroundColor = '#FFFFFF';
            _.borderColor = '#FFFFFF';
            _.barPeercentage = 0.5;
            _.titleFontSize = 16;
            _.titleColor = '#9da5cc';
            _.bodyFontSize = 16;
            _.bodyColor = '#9da5cc';
            _.hintColor = '#009900';
            _.noHintColor = '#990000';
            _.usedAllHintsColor = '#996600';
            _.gridColor = '#9da5cc';
            _.axisColor = '#212121';
            _.tickColor = '#212121';
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

