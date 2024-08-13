import m from "mithril";
import Parsley from "parsleyjs";

import Navbar from "../components/navbar";

var Login = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <div class="container">
                    <div class="block login-block">
                        <h1>Log In</h1>
                        <form class="block" id="loginForm" data-parsley-validate>
                            <input type="hidden" id="csrf_token"></input>
                            <div id="loginForm-inputs">
                                <div class="field login-field">
                                    <label class="label">Username</label>
                                    <div class="control">
                                        <input class="input" type="text" id="usernameInput" required></input>
                                    </div>
                                </div>
                                <div class="field login-field">
                                    <label class="label">Password</label>
                                    <div class="control">
                                        <input class="input" type="password" id="passwordInput" required></input>
                                    </div>
                                </div>
                            </div>
                            <div class="control">
                                <input type="submit" class="button is-link" id="signupButton" value="Log In"></input>
                            </div>
                        </form>
                        {location.href.endsWith("not_signed_in") && <p id="profile-unauthorised">Log in to view your profile!</p>}
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {

        m.request({
            method: "GET",
            url: "/get_csrf_token"
        }).then(response => {
            $("#csrf_token").val(response.csrf_token);
        });

        $("#loginForm").parsley({
            trigger: "change",
            errorsWrapper: '<div class="parsley-errors-list"></ul>',
            errorTemplate: "<p></p>"
        }).on("form:submit", () => {
            m.request({
                method: "POST",
                url: "/login",
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "username": $("#usernameInput").val(),
                    "password": $("#passwordInput").val()
                }
            }).then(response => {
                // Check if user has visited site before
                var first_visit = false;
                var URL = '#!/tutorial';
                if(localStorage.getItem('was_visited')){
                    URL = '#!/home';
                }
                first_visit = true;
                localStorage.setItem('was_visited', 1);
                console.log(first_visit);
                alert(response.msg);
                localStorage.setItem('jwt', response.access_token); // log in the user
                location.href = URL // redirect to homepage
            }).catch(e => {
                if (e.code === 401) { // unauthorized - probably invalid username/password
                    alert(e.response.msg);
                }
            });
        });
    }
}

export default Login;