import m from "mithril";

import Navbar from "../components/navbar";
import Settingsbar from "../components/settingsbar";

var settingsProfile = {
    view: () => (
        <div class="content">
            <Navbar/>
            <div class="container-wrapper">
                <div class="container">
                    <div class="block settings-block">
                        <Settingsbar/>
                        <h2>Edit Profile</h2>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            // TODO: add edit profile stuff
        }
        else{
            location.href = "#!/login?message=not_signed_in";
        }
    }
}

export default settingsProfile;