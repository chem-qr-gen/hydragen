import m from "mithril";

import Navbar from "../components/navbar";

var Landing = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="hero is-large" id="landing-splash">
                <div class="hero-body">
                    <p class="title">Mass Spectrometry Practice Made Easy</p>
                    <div class="block">
                        <a class="button is-primary mx-1" href="#!/signup">Sign Up</a>
                        <a class="button mx-1" href="#!/login">Log In</a>
                        <a class="button mx-1" href="#!/home">Try Now</a>
                    </div>
                </div>
            </div>
        </div>
    )
}

export default Landing;