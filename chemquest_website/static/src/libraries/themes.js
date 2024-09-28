function getChartStylesFromTheme(_, theme) {
    switch (theme) {
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
            _.labelFontFamily = 'Figtree'
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
            _.labelFontFamily = 'Figtree'
            _.labelFontSize = '16'
            _.labelFontWeight = '600'
            break;
    }
    return _;
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
    getChartStylesFromTheme,
    smiDrawerTheme,
}