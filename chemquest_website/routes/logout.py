from flask_login import logout_user

from chemquest_website import app

@app.route('/logout')
def logout():
    '''API for logout.'''
    logout_user()
    return {"msg": "Logout successful"}