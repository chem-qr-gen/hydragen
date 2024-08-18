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
                        <table id="leaderboard-table">
                            <thead>
                                <tr>
                                    <th scope="col">Rank</th>
                                    <th scope="col">Username</th>
                                    <th scope="col">Rating</th>
                                </tr>
                            </thead>
                            <tbody id="leaderboard-tbody"></tbody>
                        </table>
                    </div>
                </div>
            </div>
         </div>
    ),

    oncreate: () => {
        m.request({
            method: "GET",
            url: "/get_leaderboard"
        }).then(response => {
            // add in the leaderboard entries based on decreasing elo
            let users_by_rating = response;
            users_by_rating = users_by_rating.sort((a, b) => b.elo - a.elo);

            for (let id = 0; id < users_by_rating.length; id++){
                let user = users_by_rating[id];
                $("#leaderboard-tbody").append(
                    $("<tr>")
                        .append($("<th>").attr("scope", "row").text(id + 1))
                        .append($("<td>").text(user.username))
                        .append($("<td>").text(Math.round((user.elo + Number.EPSILON) * 100) / 100)) // round to 2 digits (without trailing zeroes)
                );
            }
        })
    }
}

export default Leaderboard;