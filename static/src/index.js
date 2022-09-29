import m from "mithril";
import Home from "./views/home";
import Login from "./views/login";
import Signup from "./views/signup";

var root = document.body;

m.route(root, '/home', {
    '/home': Home,
    '/login': Login,
    '/signup': Signup
})