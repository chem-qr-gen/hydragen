function checkTheme(themeName) {
    localStorage.getItem('theme', themeName)
}

function setTheme() {
	const cachedTheme = localStorage.getItem('theme');  //attempt to retrieve theme from cache
	if (!(cachedTheme)) {
		cachedTheme = 'auto';
	}
	document.documentElement.dataset['theme'] = cachedTheme;
	console.log(cachedTheme);
	return cachedTheme;
}
// Applies current theme, make sure correct button is checked, listen for user input
function applyTheme() {
	const cachedTheme = setTheme(); //Check for and apply any chached theme
	const initialTheme = cachedTheme ?? 'auto'; //Retrieve theme currently in use
	const themePicker = document.getElementById('theme-input'); 
	const a=themePicker.querySelector('input[checked]'); //Find any theme button that is currently checked
	if (!(a===null)) {
		a.removeAttribute('checked'); //Uncheck the theme button
	}
	themePicker.querySelector(`input[value="${initialTheme}"]`).setAttribute('checked', '');  //Check the theme button that matches with the current theme
	themePicker.addEventListener('change', (e) => {			//Listen for user clicking on the theme button
		const theme = e.target.value;
		document.documentElement.dataset['theme'] = theme;
		localStorage.setItem('theme', theme);
	})
}
// Runs before website loads to ensure correct theme is in use
setTheme();