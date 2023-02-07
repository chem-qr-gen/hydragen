from flask import render_template

from chemquest_website import app

@app.route('/')
def index():
    '''Returns a skeleton index page for Flask. For the actual HTML, refer to files in static/src/views.'''
    return render_template("index.html")
