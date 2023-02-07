from flask import request
from werkzeug.security import check_password_hash
from flask_jwt_extended import create_access_token

from chemquest_website import app, db

@app.route('/login', methods = ['POST'])
def login():
    '''API for login.'''
    user = db.users.find_one({"username": request.json["username"]})
    if user and check_password_hash(user["password"], request.json["password"]):
        access_token = create_access_token(user["username"])
        return {"msg": "Login successful", "access_token": access_token}
    return {"msg": "Incorrect username or password"}, 401