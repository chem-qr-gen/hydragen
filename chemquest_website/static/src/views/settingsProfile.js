import m from "mithril";

import Navbar from "../components/navbar";
import Settingsbar from "../components/settingsbar";
import CountryInput from "../components/countryInput";

var settingsProfile = {
    view: () => (
        <div class="content">
            <Navbar/>
            <div class="container-wrapper">
                <div class="container">
                    <div class="block settings-block">
                        <Settingsbar/>
                        <h2>Edit Profile</h2>
                        <form className="block" id="settingsProfileForm">
                            <input type="hidden" id="csrf_token"></input>
                            <div className="field" id="settingsProfileForm-inputs">
                                <div className="field pt-0">
                                    <label className="label">Gender</label>
                                    <div className="control">
                                        <div className="select">
                                            <select id="genderInput" required>
                                                <option disabled selected>---</option>
                                                <option>Male</option>
                                                <option>Female</option>
                                                <option>Other</option>
                                                <option>Prefer not to say</option>
                                            </select>
                                        </div>
                                    </div>
                                </div>
                                <div className="field pt-0">
                                    <label className="label">Country</label>
                                    <div className="control">
                                        <div className="select">
                                            <CountryInput id="countryInput"/>
                                        </div>
                                    </div>
                                </div>
                                <div className="field pt-0">
                                    <label className="label">Region</label>
                                    <div className="control">
                                        <input className="input" id="regionInput" type="text"/>
                                    </div>
                                </div>
                            </div>
                            <div className="control">
                                <input type="submit" className="button is-link" id="editProfileButton"
                                       value="Submit"></input>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        console.log("Test");
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            m.request({ // get csrf token
                method: "GET",
                url: "/get_csrf_token"
            }).then(response => {
                $("#csrf_token").val(response.csrf_token);
            });

            m.request({ // get user data
                method: "GET",
                url: "/get_profile",
                headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")}
            }).then(response => { // display data
                // gender select
                let genderInput = $("#genderInput")[0];
                // genderInput.selectedIndex = 1;
                for (let i = 0; i < genderInput.options.length; i++){
                    if (genderInput.options[i].text === response.gender) {
                        genderInput.selectedIndex = i;
                        break;
                    }
                }
                // country select
                let countryInput = $("#countryInput")[0];
                for (let i = 0; i < countryInput.options.length; i++){
                    if (countryInput.options[i].text === response.country) {
                        countryInput.selectedIndex = i;
                        break;
                    }
                }
                // region input
                $("#regionInput")[0].value = response.region;

                // submit form
                $("#settingsProfileForm").parsley({
                    trigger: "change",
                    errorsWrapper: '<div class="parsley-errors-list"></ul>',
                    errorTemplate: "<p></p>"
                }).on("form:submit", () => {
                    m.request({
                        method: "POST",
                        url: "/edit_profile",
                        headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")},
                        body: {
                            "_csrf_token": $("#csrf_token").val(),
                            "gender": $("#genderInput option:selected").text(),
                            "country": $("#countryInput option:selected").text(),
                            "region": $("#regionInput").val()
                        }
                    }).then(response => {
                        alert(response.msg)
                        location.reload()
                    })
                });
            });
        }
        else { // no user logged in, redirect to login page
            location.href = "#!/login?message=not_signed_in";
        }
    }
}

export default settingsProfile;
