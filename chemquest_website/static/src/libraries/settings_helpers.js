function checkTheme(themeName) {
    localStorage.getItem('theme', themeName)
}

function setTheme() {
	const cachedTheme = localStorage.getItem('theme');  //attempt to retrieve theme from cache
	if (cachedTheme) {
		document.documentElement.dataset['theme'] = cachedTheme; //set theme to the cached theme
	}
	return cachedTheme;
}

function applyTheme() {
	const cachedTheme = setTheme();
	const initialTheme = cachedTheme ?? 'auto';
	const themePicker = document.getElementById('theme-input');
	const a=themePicker.querySelector('input[checked]');
	if (!(a===null)) {
		a.removeAttribute('checked');
	}
	themePicker.querySelector(`input[value="${initialTheme}"]`).setAttribute('checked', ''); 
	themePicker.addEventListener('change', (e) => {
		const theme = e.target.value;
		document.documentElement.dataset['theme'] = theme;
		localStorage.setItem('theme', theme);
	})
}
setTheme();