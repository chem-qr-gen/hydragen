function checkTheme(themeName) {
    localStorage.getItem('theme', themeName)
}

function setTheme(themeName) {
    if (theme === 'auto') {
        delete document.documentElement.dataset['theme'];
        localStorage.removeItem('theme');
      } else {
        document.documentElement.dataset['theme'] = theme;
        localStorage.setItem('theme', theme);
      }
}

