import m from "mithril";
import Landing from "./views/landing";
import Home from "./views/home";
import Login from "./views/login";
import Signup from "./views/signup";
import Legacy from "./views/legacy";

var root = document.body;

m.route(root, '/landing', {
    '/landing': Landing,
    '/legacy': Legacy,
    '/login': Login,
    '/signup': Signup,
    '/home': Home,
})