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

(function () {
    const cachedTheme = localStorage.getItem('theme')
    if (cachedTheme) {
        document.documentElement.dataset['theme'] = cachedTheme;
    }
    const initialTheme = cachedTheme ?? 'light';
    settings_helpers.querySelector('input[checked]').removeAttribute('checked');
    settings_helpers.querySelector(`input[value="${initialTheme}"]`).setAttribute('checked', '');
    settings_helpers.addEventListener('change', (e) => {
        const theme = e.target.value;
        if (theme === 'auto') {
          delete document.documentElement.dataset['theme'];
          localStorage.removeItem('theme');
        } else {
          document.documentElement.dataset['theme'] = theme;
          localStorage.setItem('theme', theme);
        }
    })
} )();