import m from "mithril";

import Navbar from "../components/navbar";

var Landing = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="hero is-large" id="landing-splash">
                <div class="hero-body">
                    <p class="title" id="landing-title">Mass Spectrometry Practice Made Easy</p>
                    <div class="block">
                        <a class="button is-primary mx-1" id="tryNowButton" href="#!/home">Try Now</a>
                        <a class="button mx-1" id="loginButton" href="#!/login">Log In</a>
                        <a class="button mx-1" id="signupButton" href="#!/signup">Sign Up</a>
                    </div>
                </div>
                <a id="bg-src-link" href="https://www.vecteezy.com/free-videos/science">Background Source: Science Stock Videos by Vecteezy</a>
            </div>
        </div>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") !== null) { // if there is a user logged in
            $("#loginButton, #signupButton").hide();
        }
    }
}

export default Landing;