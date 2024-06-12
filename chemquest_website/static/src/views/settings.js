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
                            <div class="settings-label">
                            <h4>Change Theme</h4>
                            </div>
                            <div id="theme-input">
                                <div class="settings-control">
                                    <input type="radio" id="theme1" name="theme-selector" value="auto" onclick="setTheme(auto)" >
                                        <img id="theme-preview1"></img>
                                    </input>
                                </div>
                                <div class="settings-control">
                                    <input type="radio" id="theme2" name="theme-selector" value="light" onclick="setTheme(light)" checked >
                                        <img id="theme-preview2"></img>
                                    </input>
                                </div>
                                <div class="settings-control">
                                    <input type="radio" id="theme3" name="theme-selector" value="dark" onclick="setTheme(dark)">
                                        <img id="theme-preview3"></img>
                                    </input>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    ),
    oncreate: () => {
        (function () {
            applyTheme();
        } )();
    }
}

export default Settings;


