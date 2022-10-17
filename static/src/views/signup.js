import m from "mithril";
import $ from "cash-dom";
import validator from "email-validator";
import Pristine from "pristine";

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
                            <p class="help is-danger" id="usernameHelp"></p>
                        </div>
                    </div>
                    <div class="field">
                        <label class="label">Email</label>
                        <div class="control">
                            <input class="input" type="email" id="emailInput"></input>
                            <p class="help is-danger" id="emailHelp"></p>
                        </div>
                    </div>
                    <div class="field">
                        <label class="label">Password</label>
                        <div class="control">
                            <input class="input" type="password" id="passwordInput"></input>
                        </div>
                        <p class="help is-danger" id="passwordHelp"></p>
                    </div>
                    <div class="field">
                        <label class="label">Confirm Password</label>
                        <div class="control">
                            <input class="input" type="password" id="confirmPasswordInput"></input>
                        </div>
                        <p class="help is-danger" id="confirmPasswordHelp"></p>
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

        $("#usernameInput").on("blur change", () => {
            if ($("#usernameInput").val().length == 0) {
                $("#usernameInput").addClass("is-danger");
                $("#usernameHelp").text("Please enter a username.");
            } else {
                $("#usernameInput").removeClass("is-danger");
                $("#usernameHelp").text("");
            }
        });

        $("#emailInput").on("blur change", () => {
            if ($("#emailInput").val().length == 0) {
                $("#emailInput").addClass("is-danger");
                $("#emailHelp").text("Please enter an email.");
            } else if (!validator.validate($("#emailInput").val())) {
                $("#emailInput").addClass("is-danger");
                $("#emailHelp").text("Please enter a valid email address.");
            } else {
                $("#emailInput").removeClass("is-danger");
                $("#emailHelp").text("");
            }
        });

        $("#passwordInput").on("blur change", () => {            
            if ($("#passwordInput").val().length == 0) {
                $("#passwordInput").addClass("is-danger");
                $("#passwordHelp").text("Please enter a password.");
            } else if ($("#passwordInput").val().length < 8) {
                $("#passwordInput").addClass("is-danger");
                $("#passwordHelp").text("Minimum 8 characters.");
            } else {
                $("#passwordInput").removeClass("is-danger");
                $("#passwordHelp").text("");
            }
        });

        $("confirmPasswordInput").on("blur change", () => {
            if ($("#confirmPasswordInput").val().length == 0) {
                $("#confirmPasswordInput").addClass("is-danger");
                $("#confirmPasswordHelp").text("Please enter a password.");
            } else if ($("#confirmPasswordInput").val() !== $("#passwordInput").val()) {
                $("#confirmPasswordInput").addClass("is-danger");
                $("#confirmPasswordHelp").text("Passwords do not match.");
            } else {
                $("#confirmPasswordInput").removeClass("is-danger");
                $("#confirmPasswordHelp").text("");
            }
        });

        $("#signupButton").on("click", () => {
            $("#usernameInput").blur();
            $("#emailInput").blur();
            $("#passwordInput").blur();
            $("#confirmPasswordInput").blur();

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
                location.href = "#!/login" // redirect after successful signup
            }).catch(e => {
                if (e.code === 401) { // unauthorized - probably user already exists
                    alert(e.response.msg);
                }
            })
        });
    }
}

export default Signup;