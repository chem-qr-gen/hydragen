function getTheme() {
    var theme = document.documentElement.getAttribute('data-theme');
    if (theme == 'auto') {
        if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
            theme = 'dark';
        }
        else {
            theme = 'light';
        }
    }
    return theme;
}

export {
    getTheme,
}