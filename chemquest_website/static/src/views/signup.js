import m from "mithril";
import Parsley from "parsleyjs";

import Navbar from "../components/navbar";
import CountryInput from "../components/countryInput";

var Signup = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <div class="container" id="signup-container">
                    <div class="block">
                        <h1>Sign Up</h1>
                    </div>
                    <form class="block" id="signupForm" data-parsley-validate>
                        <input type="hidden" id="csrf_token"></input>
                        <div class="columns is-mobile is-multiline" id="signupForm-inputs">
                            <div class="field column is-half pt-0">
                                <label class="label">Username</label>
                                <div class="control">
                                    <input class="input" type="text" id="usernameInput" required></input>
                                </div>
                            </div>
                            <div class="field column is-half pt-0">
                                <label class="label">Email</label>
                                <div class="control">
                                    <input class="input" type="email" id="emailInput" required></input>
                                </div>
                            </div>
                            <div class="field column is-half pt-0">
                                <label class="label">Password</label>
                                <div class="control">
                                    <input class="input" type="password" id="passwordInput" required minlength="8"></input>
                                </div>
                            </div>
                            <div class="field column is-half pt-0">
                                <label class="label">Confirm Password</label>
                                <div class="control">
                                    <input class="input" type="password" id="confirmPasswordInput" required minlength="8" data-parsley-equalto="#passwordInput"></input>
                                </div>
                            </div>
                            <div class="field column is-one-third pt-0">
                                <label class="label">Gender</label>
                                <div class="control">
                                    <div class="select">
                                        <select id="genderInput" required>
                                            <option>---</option>
                                            <option>Male</option>
                                            <option>Female</option>
                                            <option>Other</option>
                                            <option>Prefer not to say</option>
                                        </select>
                                    </div>
                                </div>
                            </div>
                            <div class="field column is-one-third pt-0">
                                <label class="label">Country</label>
                                <div class="control">
                                    <div class="select">
                                        <CountryInput id="countryInput" />
                                    </div>
                                </div>
                            </div>
                            <div class="field column is-one-third pt-0">
                                <label class="label">Region</label>
                                <div class="control">
                                    <input class="input" type="text" id="regionInput" required></input>
                                </div>
                            </div>
                        </div>
                        <div class="control">
                            <input type="submit" class="button is-link" id="signupButton" value="Sign Up"></input>
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

        // TODO: make the form dynamically check if username and email already exist (should be possible with parsley)
        $("#signupForm").parsley({
            trigger: "change",
            errorsWrapper: '<div class="parsley-errors-list"></ul>',
            errorTemplate: "<p></p>"
        }).on("form:submit", () => {
            m.request({
                method: "POST",
                url: "/signup",
                body: {
                    "_csrf_token": $("#csrf_token").val(),
                    "username": $("#usernameInput").val(),
                    "email": $("#emailInput").val(),
                    "password": $("#passwordInput").val(),
                    "gender": $("#genderInput option:selected").text(),
                    "country": $("#countryInput option:selected").text(),
                    "region": $("#regionInput").val()
                }
            }).then(response => {
                alert(response.msg);
                location.href = "#!/login" // redirect after successful signup
            }).catch(e => {
                if (e.code === 401) { // unauthorized - probably user already exists
                    alert(e.response.msg);
                }
            });
            return false;
        });
    }
}

export default Signup;