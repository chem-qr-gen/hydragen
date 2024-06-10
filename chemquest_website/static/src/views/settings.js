import m from "mithril";
import Parsley from "parsleyjs";

import Navbar from "../components/navbar";

var Settings = {
    view: () => (
        <div class="content">
            <Navbar/>
            <div class="container-wrapper">
                <div class="container">
                    <div class="block settings-block">
                        <h2>Graphics</h2>
                        <div class="field settings-field">
                            <div class="settings-label">
                                <h4>Play Animations</h4>
                            </div>
                            <div class="settings-input">
                                <input type="checkbox"></input>
                                <span class="toggle-switch"></span>
                            </div>
                        </div>
                        <div class="field settings-field">
                            <h4>Change Theme</h4>
                            <div class="settings-input">
                                <div class="settings-control">
                                    <input type="radio" name="theme-selector" value="0"></input>
                                    <span id="theme1"></span>
                                </div>
                                <div class="settings-control">
                                    <input type="radio" name="theme-selector" value="1"></input>
                                    <span id="theme2"></span>
                                </div>
                                <div class="settings-control">
                                    <input type="radio" name="theme-selector" value="2"></input>
                                    <span id="theme3"></span>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    )
}

export default Settings;