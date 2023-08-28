from flask import request
from werkzeug.security import generate_password_hash
from flask_mail import Message

from chemquest_website import app, engine, meta, mail

@app.route('/signup', methods = ['POST'])
def signup():
    '''API for signup.'''

    users_table = meta.tables["users"]

    # check if user already exists
    with engine.connect() as conn:
        existing_user = conn.execute(
            users_table.select().where(users_table.c.username == request.json["username"])
        ).fetchone()

    # existing_user = db.users.find_one({"username": request.json["username"]}) # check if a user already exists
    if existing_user:
        return {"msg": "Username already exists"}, 401
    
    # user doesn't already exist, so add them to the database
    new_user = {
        "username": request.json["username"],
        "email": request.json["email"],
        "password": generate_password_hash(request.json["password"]),
        "gender": request.json["gender"],
        "country": request.json["country"],
        "region": request.json["region"],
        "elo": 1000.0
    }
    with engine.connect() as conn:
        conn.execute(users_table.insert().values(**new_user))

    # db.users.insert_one(new_user).inserted_id

    # send confirmation email to user
    email_msg = Message(
        subject="ChemQuest - Signup Successful", 
        recipients=[new_user["email"]], 
        body=f"Hello {new_user['username']},\n\nCongratulations! You have successfully signed up for ChemQuest."
    )
    mail.send(email_msg)
    return {"msg": "Signup successful"}
