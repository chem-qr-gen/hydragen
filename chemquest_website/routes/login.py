from box import Box
from flask import request
from flask_login import login_user
from werkzeug.security import check_password_hash

from chemquest_website import app, engine, meta, User

@app.route('/login', methods = ['POST'])
def login():
    '''API for login. If login is unsuccessful, 401 will be returned.'''

    # retrieve the user from the database
    users_table = meta.tables["users"]
    with engine.connect() as conn:
        user = conn.execute(
            users_table.select().where(users_table.c.username == request.json["username"])
        ).fetchone()

    # check the password. if it's correct, login the user using flask-login. otherwise, return an error.
    if user:
        user_dict = Box(user._mapping)
        if check_password_hash(user_dict["password"], request.json["password"]):
            user_obj = User(user_dict.username)
            login_user(user_obj)
            return {"msg": "Login successful"}
        
    return {"msg": "Incorrect username or password"}, 401