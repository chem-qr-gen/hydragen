import m from "mithril";

import Navbar from "../components/navbar";


var Landing = {
    view: () => (
        <div class="content">
            <Navbar />
            <video autoplay muted loop plays-inline class="background-video">
                <source src="../static/images/landing.mp4" type="video/mp4"></source>
            </video>
            <div class="unbounded-wrapper">
                <div class="hero is-large" id="landing-splash">
                    <div class="hero-body" id="title-container">
                        <p class="title" id="landing-title">Mass Spectrometry Practice Made Easy</p>
                        <div class="block">
                            <a class="button is-primary mx-1" id="tryNowButton" href="#!/tutorial">Try Now</a>
                            <a class="button mx-1" id="loginButton" href="#!/login">Log In</a>
                            <a class="button mx-1" id="signupButton" href="#!/signup">Sign Up</a>
                        </div>
                    </div>
                    <div>
                        <div id="statistics-div">
                            <h4>Statistics</h4>
                            <p>
                                <span className="statistics-label">Users: </span>
                                <span id="statistics-numUsersSpan"></span>
                            </p>
                            <p>
                                <span className="statistics-label">Attempts: </span>
                                <span id="statistics-numAttemptsSpan"></span>
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        // attempt to request /current_user to check if user is logged in
        // if user is logged in, hide login and signup buttons
        m.request({
            method: "GET",
            url: "/current_user"
        }).then(response => {
            if (response.username) {
                $("#loginButton, #signupButton").hide();
            }
        }).catch(() => {
            console.log("No user logged in");
        });

        m.request({ // get site statistics
            method: "GET",
            url: "/get_site_statistics"
        }).then(response => { // display data
            $("#statistics-numUsersSpan").text(response.users_count);
            $("#statistics-numAttemptsSpan").text(response.attempts_count);
        })
        $("#tryNowButton").on("Click", async () => {
            localStorage.setItem("was_visited", 1);
        }) 
    }
}

export default Landing;