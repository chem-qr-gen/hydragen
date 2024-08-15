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
                        <h2>Edit Profile (WIP)</h2>
                        <form className="block" id="settingsProfileForm">
                            <div id="settingsProfileForm-inputs">
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
                        </form>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            m.request({ // get user data
                method: "GET",
                url: "/get_profile",
                headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")}
            }).then(response => { // display data
                // gender select
                /*let genderInput = $("#genderInput");
                for (let i = 0; i < genderInput.options.length; i++){
                    if (genderInput.options[i].text === response.gender) {
                        genderInput.selectedIndex = i;
                        break;
                    }
                }
                // country select
                let countryInput = $("#countryInput");
                for (let i = 0; i < countryInput.options.length; i++){
                    if (countryInput.options[i].text === response.country) {
                        countryInput.selectedIndex = i;
                        break;
                    }
                }
                // region input
                $("#regionInput").value = response.region;
                */
            });
        }
        else { // no user logged in, redirect to login page
            location.href = "#!/login?message=not_signed_in";
        }
    }
}

export default settingsProfile;