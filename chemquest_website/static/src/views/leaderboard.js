import m from "mithril";

import Navbar from "../components/navbar";

var Leaderboard = {
    view: () => (
         <div class="content">
            <Navbar/>
            <div class="container-wrapper">
                <div class="container">
                    <div class="block leaderboard-block">
                        <h1>Leaderboard</h1>
                    </div>
                </div>
            </div>
         </div>
    ),
    oncreate: () => {

    }
}

export default Leaderboard;