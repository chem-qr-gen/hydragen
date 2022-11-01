import m from "mithril";

import Navbar from "../components/navbar";

var Landing = {
    view: () => (
        <div class="content">
            <Navbar />
            <div class="hero is-fullheight-with-navbar" id="landing-splash">
                <div class="hero-body">
                    <p class="title">Mass Spectrometry Practice Made Easy</p>
                    <div class="block">
                        <a class="button is-primary" href="#!/signup">Sign Up</a>
                        <a class="button" href="#!/login">Log In</a>
                        <a class="button" href="#!/home">Try Now</a>
                    </div>
                </div>
            </div>
            <div class="container">
                <h1>test</h1>
            </div>
        </div>
    )
}

export default Landing;