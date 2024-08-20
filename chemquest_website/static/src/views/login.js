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

        $("#loginForm").parsley({
            trigger: "change",
            errorsWrapper: '<div class="parsley-errors-list"></ul>',
            errorTemplate: "<p></p>"
        }).on("form:submit", () => {
            m.request({
                method: "POST",
                url: "/login",
                body: {
                    "username": $("#usernameInput").val(),
                    "password": $("#passwordInput").val()
                }
            }).then(response => {
                // Check if user has visited site before
                var first_visit = false;
                var redirect_url;
                if(localStorage.getItem('was_visited')){
                    redirect_url = '#!/home';
                } else {
                    redirect_url = '#!/tutorial';
                }
                first_visit = true;
                localStorage.setItem('was_visited', 1);
                console.log(first_visit);
                alert(response.msg);
                location.href = redirect_url // redirect to homepage
            }).catch(e => {
                if (e.code === 401) { // unauthorized - probably invalid username/password
                    alert(e.response.msg);
                }
            });
            return false;
        });
    }
}

export default Login;