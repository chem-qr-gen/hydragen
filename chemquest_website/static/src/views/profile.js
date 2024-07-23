import m from "mithril";

import Navbar from "../components/navbar";

var Profile = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <div class="container">
                    <div class="block profile-block">
                        <h1>Profile</h1>
                        <div class="field profile-field">
                            <p>
                                <span class="profile-label">Username: </span>
                                <span id="usernameSpan"></span>
                            </p>
                            <p>
                                <span class="profile-label">Rating: </span>
                                <span id="ratingSpan"></span>
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            // TODO: display user profile
        }
        else { // no user logged in, redirect to login page
            alert("Please log in to view your profile!");
            location.href = "#!/login";
        }
    }
}

export default Profile;