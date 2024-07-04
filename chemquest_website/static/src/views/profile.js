import m from "mithril";

import Navbar from "../components/navbar";

var Profile = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="container-wrapper">
                <div class="container">
                    <h1>Profile Placeholder</h1>
                    <p>Username: <span id="usernameSpan"></span></p>
                    <p>Rating: <span id="ratingSpan"></span></p>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        // TODO: code to retrieve user data and display it
    }
}