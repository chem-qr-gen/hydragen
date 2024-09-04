import m from "mithril";

var Settingsbar = {
    view: () => (
        <div class="settingsbar">
            <ul class="settingsbar-links-list">
                <li><a href="#!/settings">General</a></li>
                <li><a href="#!/settings/profile">Profile</a></li>
            </ul>
        </div>
    ),
    oncreate: () => {
        // if not logged in, hide settings bar
        m.request({
            method: "GET",
            url: "/current_user"
        }).catch(() => {
            $(".settingsbar").toggle();
        });
    }
};

export default Settingsbar;