from flask import request
from werkzeug.security import check_password_hash
from flask_jwt_extended import create_access_token

from chemquest_website import app, engine, meta

@app.route('/login', methods = ['POST'])
def login():
    '''API for login. If login is unsuccessful, 401 will be returned.'''

    # retrieve the user from the database
    users_table = meta.tables["users"]
    with engine.connect() as conn:
        user = conn.execute(
            users_table.select().where(users_table.c.username == request.json["username"])
        ).fetchone()
    # user = db.users.find_one({"username": request.json["username"]})

    # check the password. if it's correct, return an access token. otherwise, return an error.
    if user:
        user = user._mapping
        if check_password_hash(user["password"], request.json["password"]):
            access_token = create_access_token(user["username"])
            return {"msg": "Login successful", "access_token": access_token}
    return {"msg": "Incorrect username or password"}, 401