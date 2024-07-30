import m from "mithril";

var Settingsbar = {
    view: () => (
        <div class="settingsbar">
            <ul class="settingsbar-links-list">
                <li><a href="#!/settings">General</a></li>
                <li><a>Profile</a></li>
            </ul>
        </div>
    ),
    oncreate: () => {
        if (localStorage.getItem("jwt") === null) { // hide if not logged in
            $(".settingsbar").toggle();
        }
    }
};

export default Settingsbar;