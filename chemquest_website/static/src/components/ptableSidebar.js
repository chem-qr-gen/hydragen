import m from "mithril";

// TODO: make this responsive
var PtableSidebar = {
    view: () => (
        <div class="ptable-sidebar-wrapper ptable-sidebar-closed">
            <div class="ptable-sidebar">
                <img src="/static/images/Periodic_table_simple_en.svg" style="margin: 32px;"></img>
            </div>
            <div class="ptable-pulltab">
                ◀ Periodic Table
            </div>
        </div>
    ),
    oncreate: () => {
        $(".ptable-pulltab").on("click", () => {
            $(".ptable-sidebar-wrapper").toggleClass("ptable-sidebar-closed");
        });
    }
};

export default PtableSidebar;