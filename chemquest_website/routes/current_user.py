from flask_login import current_user, login_required

from chemquest_website import app

@app.route('/current_user')
@login_required
def logged_in_user():
    return {"username": current_user.get_id()}