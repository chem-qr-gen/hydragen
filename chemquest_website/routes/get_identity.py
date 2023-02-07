from flask_jwt_extended import get_jwt_identity, jwt_required

from chemquest_website import app

@app.route('/get_identity')
@jwt_required()
def get_identity():
    '''Returns the username of the current logged-in user, if logged in. Otherwise throws a 401.'''
    return {"identity": get_jwt_identity()}
