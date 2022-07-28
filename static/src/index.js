import m from "mithril";
import Home from "./views/home";

var root = document.body;

m.route(root, '/home', {
    '/home': Home
})