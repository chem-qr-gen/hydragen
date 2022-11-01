import m from "mithril";
import Landing from "./views/landing";
import Home from "./views/home";
import Login from "./views/login";
import Signup from "./views/signup";

var root = document.body;

m.route(root, '/landing', {
    '/landing': Landing,
    '/home': Home,
    '/login': Login,
    '/signup': Signup
})