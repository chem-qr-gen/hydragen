import m from "mithril";

import Navbar from "../components/navbar";

var Leaderboard = {
    view: () => (
         <div class="content">
            <Navbar/>
            <div class="unbounded-wrapper">
                <div class="container">
                    <div class="block leaderboard-block">
                        <h1>Leaderboard</h1>
                        <ol id="leaderboard-list"></ol>
                    </div>
                </div>
            </div>
         </div>
    ),

    oncreate: () => {
        m.request({
            method: "GET",
            url: "/get_leaderboard"
        }).then(response => { // add in the leaderboard entries based on decreasing elo
            let users_by_rating = response;
            users_by_rating = users_by_rating.sort((a, b) => b.elo - a.elo);
            let user;
            for (user of users_by_rating){
                $("#leaderboard-list").append("<li><span>" + user.username + "</span><span>" + user.elo + "</span></li>")
            }
        })
    }
}

export default Leaderboard;