from flask import request
from werkzeug.security import generate_password_hash
from flask_mail import Message

from chemquest_website import app, db, mail

@app.route('/signup', methods = ['POST'])
def signup():
    '''API for signup.'''
    existing_user = db.users.find_one({"username": request.json["username"]}) # check if a user already exists
    if existing_user:
        return {"msg": "Username already exists"}, 401
    new_user = {
        "username": request.json["username"],
        "email": request.json["email"],
        "password": generate_password_hash(request.json["password"]),
        "gender": request.json["gender"],
        "country": request.json["country"],
        "region": request.json["region"],
        "elo": 1000.0
    }
    db.users.insert_one(new_user).inserted_id
    email_msg = Message(
        subject="ChemQuest - Signup Successful", 
        recipients=[new_user["email"]], 
        body=f"Hello {new_user['username']},\n\nCongratulations! You have successfully signed up for ChemQuest."
    )
    mail.send(email_msg)
    return {"msg": "Signup successful"}
