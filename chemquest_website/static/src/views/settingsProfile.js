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
        m.request({ // checks if there is a user logged in and gets user data
            method: "GET",
            url: "/get_profile"
        }).then(response => {
            // display data

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
                    body: {
                        "gender": $("#genderInput option:selected").text(),
                        "country": $("#countryInput option:selected").text(),
                        "region": $("#regionInput").val()
                    }
                }).then(response => {
                    alert(response.msg)
                    location.reload()
                })
            });
        }).catch(() => {
            location.href = "#!/login?message=not_signed_in";
        })
    }
}

export default settingsProfile;
