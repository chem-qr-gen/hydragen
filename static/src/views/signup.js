import m from "mithril";
import $ from "cash-dom";

import Navbar from "../components/navbar";

var Signup = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container">
                <div class="block">
                    <h1>Sign Up</h1>
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
                        <label class="label">Email</label>
                        <div class="control">
                            <input class="input" type="email" id="emailInput"></input>
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
                            <button class="button is-link" id="signupButton">Sign Up</button>
                        </div>
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
        $("#signupButton").on("click", () => {
            console.log("test");
            m.request({
                method: "POST",
                url: "/signup",
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "username": $("#usernameInput").val(),
                    "email": $("#emailInput").val(),
                    "password": $("#passwordInput").val()
                }
            }).then(response => {
                alert(response.msg);
                location.href = "#!/home"
            }).catch(e => {
                if (e.code === 401) {
                    alert(e.response.msg);
                }
            })
        });
    }
}

export default Signup;