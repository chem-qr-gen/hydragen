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
            m.request({ // get user data
                method: "GET",
                url: "/get_profile",
                headers: {"Authorization": "Bearer " + localStorage.getItem("jwt")}
            }).then(response => { // display data
                $("#usernameSpan").text(response.username);
                $("#ratingSpan").text(response.elo);
            })
        }
        else { // no user logged in, redirect to login page
            location.href = "#!/login?message=not_signed_in";
        }
    }
}

export default Profile;