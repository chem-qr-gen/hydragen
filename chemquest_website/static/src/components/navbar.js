import m from "mithril";

var Navbar = {
    view: () => (
        <nav class="navbar">
            <div class="navbar-brand">
                <a class="navbar-item is-size-5" href="#!/home">
                    <img class="logo-horizontal" />
                </a>
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
        m.request({ // get the current logged in user
            method: "GET",
            url: "/current_user"
        }).then(response => {
            // hide signup and login nav items
            $("#signupLink").toggle();
            $("#loginLink").toggle();
            // show profile and logout nav items
            $("#usernameLink").toggle();

            // set the username in the navbar
            $("#usernameLink .navbar-link").text("Logged in as " + response.username);
        });
        
        // logout functionality
        $("#logoutLink").on("click", () => {
            m.request({
                method: "GET",
                url: "/logout"
            }).then(() => {
                // reload the page after logging out
                location.reload();
            })
        });
    }
};

export default Navbar;