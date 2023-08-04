import m from "mithril";
import Parsley from "parsleyjs";

import Navbar from "../components/navbar";

var Login = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <div class="container">
                    <div class="block">
                        <h1>Log In</h1>
                    </div>
                    <form class="block" id="loginForm" data-parsley-validate>
                        <input type="hidden" id="csrf_token"></input>
                        <div class="field">
                            <label class="label">Username</label>
                            <div class="control">
                                <input class="input" type="text" id="usernameInput" required></input>
                            </div>
                        </div>
                        <div class="field">
                            <label class="label">Password</label>
                            <div class="control">
                                <input class="input" type="password" id="passwordInput" required></input>
                            </div>
                        </div>
                        <div class="field">
                            <div class="control">
                                <input type="submit" class="button is-link" id="signupButton" value="Log In"></input>
                            </div>
                        </div>
                    </form>
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
                alert(response.msg);
                localStorage.setItem('jwt', response.access_token); // log in the user
                location.href = "#!/home" // redirect to homepage
            }).catch(e => {
                if (e.code === 401) { // unauthorized - probably invalid username/password
                    alert(e.response.msg);
                }
            });
        });
    }
}

export default Login;