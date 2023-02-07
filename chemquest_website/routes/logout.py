from flask import jsonify
from flask_jwt_extended import unset_jwt_cookies

from chemquest_website import app

@app.route('/logout')
def logout():
    '''API for logout.'''
    response = jsonify({"msg": "Logout successful"})
    unset_jwt_cookies(response)
    return response