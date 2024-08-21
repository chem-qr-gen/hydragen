from flask_login import login_required, current_user

from chemquest_website import app

@app.route('/get_identity')
@login_required
def get_identity():
    '''Returns the username of the current logged-in user, if logged in. Otherwise throws a 401.'''
    return {"identity": current_user.get_id()}
