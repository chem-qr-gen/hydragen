import m from "mithril";
import Landing from "./views/landing";
import Home from "./views/home";
import Login from "./views/login";
import Signup from "./views/signup";
import Tutorial from "./views/tutorial";
import Settings from "./views/settings";
import Profile from "./views/profile";

var root = document.body;

m.route(root, '/landing', {
    '/landing': Landing,
    '/login': Login,
    '/signup': Signup,
    '/home': Home,
    '/tutorial': Tutorial,
    '/settings': Settings,
    '/profile': Profile,
})
