import m from "mithril";

var Navbar = {
    view: () => (
        <nav class="navbar">
            <div class="navbar-brand">
                <a class="navbar-item is-size-5" href="#!/home"><b><span class="brand">Hydragen.io</span></b></a>
                <div class="navbar-item is-size-5" href="#!/home"><b>Hydragen.io</b></div>
                <a role="button" class="navbar-burger" aria-label="menu" aria-expanded="false" data-target="navbarContent">
                    <span aria-hidden="true"></span>
                    <span aria-hidden="true"></span>
                    <span aria-hidden="true"></span>
                </a>
            </div>
            <div id="navbarContent" class="navbar-menu">
                <div class="navbar-end">
                    <a class="navbar-item" href="#!/signup" id="signupLink">Sign Up</a>
                    <a class="navbar-item" href="#!/login" id="loginLink">Log In</a>
                    <a class="navbar-item" href="#!/settings" id="settingsLink"><img id="settings-logo"></img></a>
                    <div class="navbar-item has-dropdown is-hoverable" id="usernameLink" style="display: none;">
                        <a class="navbar-link">Logged in as </a>
                        <div class="navbar-dropdown">
                            <a class="navbar-item" href="#!/profile">Profile</a>
                            <a class="navbar-item" href="#" id="logoutLink">Log Out</a>
                        </div>
                    </div>
                </div>
            </div>
        </nav>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            m.request({ // get the username
                method: "GET",
                url: "/get_identity",
                headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")}
            }).then(response => {
                // hide signup and login nav items
                $("#signupLink").toggle();
                $("#loginLink").toggle();
                // show profile and logout nav items
                $("#usernameLink").toggle();

                // set the username in the navbar
                $("#usernameLink .navbar-link").text("Logged in as " + response.identity);
            });
        }
        
        // logout functionality
        $("#logoutLink").on("click", () => {
            m.request({
                method: "GET",
                url: "/logout",
                headers: {
                    "Authorization": "Bearer " + localStorage.getItem("jwt")
                }
            }).then(() => {
                localStorage.removeItem('jwt');
                location.reload();
            })
        });
    }
};

export default Navbar;