import m from "mithril";
import {
  AudioPlayer, 
  randomMusic,
  nextMusic,
  musicList } from "./audioPlayer";
// import {load} from "../../.pnp.loader.mjs";

var Navbar = {
  view: () => (
    <nav class="navbar">
      <div class="navbar-brand">
        <a class="navbar-item is-size-5" href="#!/home">
          <img class="logo-horizontal"/>
        </a>
        <a role="button" class="navbar-burger" aria-label="menu" aria-expanded="false" data-target="navbarContent">
          <span aria-hidden="true"></span>
          <span aria-hidden="true"></span>
          <span aria-hidden="true"></span>
        </a>
      </div>
      <div id="navbarContent" class="navbar-menu">
        <div class="navbar-end">
          <a className="navbar-item" href="#!/signup" id="signupLink">Sign Up</a>
          <a className="navbar-item" href="#!/login" id="loginLink">Log In</a>
          <a className="navbar-item" id="musicToggle">
            <span id="music_icon" className="material-symbols-outlined navbar-icon">music_note</span>
          </a>
          <a className="navbar-item" href="#!/leaderboard" id="leaderboardLink">
            <span className="material-symbols-outlined navbar-icon">leaderboard</span>
          </a>
          <a className="navbar-item" href="#!/settings" id="settingsLink">
            <span className="material-symbols-outlined navbar-icon">settings</span>
          </a>
          <div class="navbar-item has-dropdown is-hoverable" id="usernameLink" style="display: none;">
            <a className="navbar-link">Logged in as </a>
            <div class="navbar-dropdown">
              <a className="navbar-item" href="#!/profile">Profile</a>
              <a className="navbar-item" href="#" id="logoutLink">Log Out</a>
            </div>
          </div>
        </div>
      </div>
    </nav>
  ),
  oncreate: () => {
    m.request({ // get the current logged in user
      method: "GET",
      url: "/current_user"
    }).then(response => {
      // hide signup and login nav items
      $("#signupLink").toggle();
      $("#loginLink").toggle();
      // show profile and logout nav items
      $("#usernameLink").toggle();

      // set the username in the navbar
      $("#usernameLink .navbar-link").text("Logged in as " + response.username);
    });

    // logout functionality
    $("#logoutLink").on("click", () => {
      m.request({
        method: "GET",
        url: "/logout"
      }).then(() => {
        // reload the page after logging out
        location.reload();
      })
    });
    var music = new Object();
    document.addEventListener('DOMContentLoaded', () => {
      // if (localStorage.getItem("music_play") === null) {
      localStorage.setItem("music_play", "0");
      // }
      music = randomMusic();
      console.log(music);
      var currentMusic = "get-audio/" + music.playMusic;
      AudioPlayer.load(currentMusic);
      music = nextMusic(music.nextInstrument);
      console.log(music);
      // if (localStorage.getItem("music_play") === "1") {
      //   $("#music_icon")[0].textContent = "music_note";
      //   AudioPlayer.play();
      // } else {
      $("#music_icon")[0].textContent = "music_off";
      //   AudioPlayer.pause();
      // }
      console.log("Music Time: " + AudioPlayer.audio.currentTime);
    });

    AudioPlayer.audio.addEventListener("ended", () => {
      console.log("music ended");
      AudioPlayer.audio.currentTime = 0;
      var currentMusic = "get-audio/" + music.playMusic;
      AudioPlayer.load(currentMusic);
      AudioPlayer.play();
      music = nextMusic(music.nextInstrument);
    });

    $("#music_icon").on("click", () => {
      if (localStorage.getItem("music_play") === "1") {
        localStorage.setItem("music_play", "0");
        $("#music_icon")[0].textContent = "music_off";
        AudioPlayer.pause();
        console.log("Music Time: " + AudioPlayer.audio.duration);
      } else {
        localStorage.setItem("music_play", "1");
        $("#music_icon")[0].textContent = "music_note";
        AudioPlayer.play();
        console.log("Music Time: " + AudioPlayer.audio.duration);
      }
    });
  }
};

export default Navbar;
