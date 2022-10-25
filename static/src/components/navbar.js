import m from "mithril";

var Navbar = {
    view: () => (
        <nav class="navbar">
            <div class="navbar-brand">
                <a class="navbar-item is-size-5" href="#!/home"><b>ChemQuest</b></a>
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
                $("#signupLink").toggle();
                $("#loginLink").text(response.identity + " (Log out)");
                $("#loginLink").attr("href", "#");
                $("#loginLink").on("click", () => { // turn the login link into a logout link
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
            });
        }
        
    }
};

export default Navbar;