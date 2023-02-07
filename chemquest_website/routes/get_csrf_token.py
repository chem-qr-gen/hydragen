from chemquest_website import app, csrf

@app.route('/get_csrf_token')
def get_csrf_token():
    '''CSRF token access point for CSRF protection.'''
    return {'csrf_token': csrf._get_token()}