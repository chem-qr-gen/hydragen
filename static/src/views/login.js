import m from "mithril";
import $ from "cash-dom";

var Login = {
    view: () => (
        <div class="content container">
            <div class="block">
                <h1>Log In</h1>
            </div>
            <div class="block">
                <input type="hidden" id="csrf_token"></input>
                <div class="field">
                    <label class="label">Username</label>
                    <div class="control">
                        <input class="input" type="text" id="usernameInput"></input>
                    </div>
                </div>
                <div class="field">
                    <label class="label">Password</label>
                    <div class="control">
                        <input class="input" type="password" id="passwordInput"></input>
                    </div>
                </div>
                <div class="field">
                    <div class="control">
                        <button class="button is-link" id="loginButton">Log In</button>
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
        $("#loginButton").on("click", () => {
            console.log("test");
            m.request({
                method: "POST",
                url: "/login",
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "username": $("#usernameInput").val(),
                    "password": $("#passwordInput").val()
                }
            }).then(response => {
                alert(response.msg)
            }).catch(e => {
                alert(e.response.msg)
            })
        });
    }
}

export default Login;