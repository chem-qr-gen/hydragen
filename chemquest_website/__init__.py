import datetime, json, random
from werkzeug.security import check_password_hash, generate_password_hash
from flask import Flask, redirect, render_template, request, url_for, jsonify
from flask_jwt_extended import JWTManager, create_access_token, get_jwt, get_jwt_identity, jwt_required, unset_jwt_cookies
from flask_mail import Mail, Message
from flask_seasurf import SeaSurf
from pymongo import MongoClient

import chemquest_website.elo as elo
#import mcq_generator as mcq
from chemquest_website.chemquest_secrets import * # secret addresses, change or remove for local testing

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
app.config['SECRET_KEY'] = csrf_secret
app.config['JWT_SECRET_KEY'] = jwt_secret
app.config['MAIL_SERVER'] = mail_server_secret
app.config['MAIL_DEFAULT_SENDER'] = mail_sender_secret
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USERNAME'] = mail_username_secret
app.config['MAIL_PASSWORD'] = mail_password_secret
#app.config['JWT_ACCESS_TOKEN_EXPIRES'] = datetime.timedelta(days = 1)
jwt = JWTManager(app)
mail = Mail(app)
csrf = SeaSurf(app)
client = MongoClient(mongo_address)
db = client.chemquest_db

#molecules_df = mcq.init()
#fpe = mcq.start_engine()

@app.route('/get_csrf_token')
def get_csrf_token():
    '''CSRF token access point for CSRF protection.'''
    return {'csrf_token': csrf._get_token()}

@app.route('/')
def index():
    '''Returns a skeleton index page for Flask. For the actual HTML, refer to files in static/src/views.'''
    return render_template("index.html")

@app.route('/ms_questions_new')
def ms_questions_new():
    '''Returns a set of MS data.'''
    question_id = request.args.get('id')
    if question_id == "random":
        return redirect(url_for('ms_questions_new', id = random.randint(0, db.ms_data.count_documents({}) - 1)))
    question_response = db.ms_data.find_one({"qid": int(question_id)})
    question_response.pop("_id")
    return question_response

#@app.route('/generate_mcq')
#def generate_mcq():
#    '''Returns a set of similar molecules to the input SMILES for use in MCQ questions.'''
#    input_smiles = request.args.get('input_smiles')
#    return mcq.get_similar_compounds(molecules_df, fpe, input_smiles)

@app.route('/record_attempt', methods = ['POST'])
@jwt_required(optional = True)
def record_attempt():
    '''Record a question attempt made by the user, and relevant stats.'''
    question_id = request.json['id']
    answer = request.json['answer']
    is_correct = request.json['isCorrect']

    jwt_identity = get_jwt_identity()
    if jwt_identity:
        user_dict = db.users.find_one({"username": jwt_identity})
        question_dict = db.ms_data.find_one({"qid": question_id})

        # TODO: implement elo system based on partial scores with wrong MCQ answers. Elo changes will be calculated when "next question" clicked or navigated away.
        # CURRENT: Elo is calculated after every attempt, treating it as a full win (1) or loss (0).
        player_old_elo = user_dict["elo"] if "elo" in user_dict else 1000.0
        player_new_elo = elo.rate_single(player_old_elo, question_dict["difficulty"], float(is_correct), 1.0 - float(is_correct))[0]

        this_attempt = {
            "timestamp": datetime.datetime.utcnow(),
            "question_id": question_id,
            "answer": answer,
            "is_correct": is_correct,
            "player_old_elo": player_old_elo,
            "player_new_elo": player_new_elo,
        }
        
        new_attempts = user_dict["attempts"] if "attempts" in user_dict else []
        new_attempts.append(this_attempt)
        
        db.users.update_one(
            {"_id": user_dict["_id"]},
            {"$set": {
                "attempts": new_attempts,
                "elo": player_new_elo
            }},
            upsert = False)
    
        return this_attempt
    return {"msg": "Not logged in, response not recorded."}


@app.route('/login', methods = ['POST'])
def login():
    '''API for login.'''
    user = db.users.find_one({"username": request.json["username"]})
    if user and check_password_hash(user["password"], request.json["password"]):
        access_token = create_access_token(user["username"])
        return {"msg": "Login successful", "access_token": access_token}
    return {"msg": "Incorrect username or password"}, 401

@app.route('/logout')
def logout():
    '''API for logout.'''
    response = jsonify({"msg": "Logout successful"})
    unset_jwt_cookies(response)
    return response

@app.route('/signup', methods = ['POST'])
def signup():
    '''API for signup.'''
    existing_user = db.users.find_one({"username": request.json["username"]})
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

@app.route('/get_identity')
@jwt_required()
def get_identity():
    '''Returns the username of the current logged-in user, if logged in. Otherwise throws a 401.'''
    return {"identity": get_jwt_identity()}


if __name__ == "__main__":
    app.run(debug = True, host = "0.0.0.0")
