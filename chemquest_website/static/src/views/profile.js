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
                                <span className="profile-label">Username: </span>
                                <span id="profile-usernameSpan"></span>
                            </p>
                            <p>
                                <span className="profile-label">Rating: </span>
                                <span id="profile-ratingSpan"></span>
                            </p>
                            <p>
                                <span className="profile-label">Email: </span>
                                <span id="profile-emailSpan"></span>
                            </p>
                            <p>
                                <span className="profile-label">Gender: </span>
                                <span id="profile-genderSpan"></span>
                            </p>
                            <p>
                                <span className="profile-label">Country: </span>
                                <span id="profile-countrySpan"></span>
                            </p>
                            <p>
                                <span className="profile-label">Region: </span>
                                <span id="profile-regionSpan"></span>
                            </p>
                        </div>
                        <a class="button is-link" href="#!/settings/profile">Edit Profile</a>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        // this will throw a 403 (forbidden) if the user is not logged in
        m.request({
            method: "GET",
            url: "/get_profile"
        }).then(response => {
            $("#profile-usernameSpan").text(response.username);
            $("#profile-ratingSpan").text(response.elo);
            $("#profile-emailSpan").text(response.email);
            $("#profile-genderSpan").text(response.gender);
            $("#profile-countrySpan").text(response.country);
            $("#profile-regionSpan").text(response.region);
        }).catch(() => { // if the user is not logged in, redirect to login page
            location.href = "#!/login?message=not_signed_in";
        });
    }
}

export default Profile;