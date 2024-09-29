import m from "mithril";
import Landing from "./views/landing";
import Home from "./views/home";
import Login from "./views/login";
import Signup from "./views/signup";
import Tutorial from "./views/tutorial";
import Settings from "./views/settings";
import Profile from "./views/profile";
import settingsProfile from "./views/settingsProfile";
import Leaderboard from "./views/leaderboard";
import {AudioPlayer} from "./components/audioPlayer";

var root = document.body;

// m.route(root, '/landing', {
//     '/landing': Landing,
//     '/login': Login,
//     '/signup': Signup,
//     '/home': Home,
//     '/tutorial': Tutorial,
//     '/leaderboard': Leaderboard,
//     '/settings': Settings,
//     '/settings/profile': settingsProfile,
//     '/profile': Profile,
// })
m('div', [
  m(AudioPlayer),  // The global audio player
  m('div', {id: 'content'}, [
    m.route(root, '/landing', {
      '/landing': Landing,
      '/login': Login,
      '/signup': Signup,
      '/home': Home,
      '/tutorial': Tutorial,
      '/leaderboard': Leaderboard,
      '/settings': Settings,
      '/settings/profile': settingsProfile,
      '/profile': Profile,
    })
  ])
]);
