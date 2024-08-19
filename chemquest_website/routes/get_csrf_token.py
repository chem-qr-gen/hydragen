from flask import abort

from chemquest_website import app

@app.route('/get_csrf_token')
def get_csrf_token():
    '''CSRF token access point for CSRF protection.'''
    # return 404, deprecated functionality
    abort(404)