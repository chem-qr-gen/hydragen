import datetime
from flask import request
from flask_jwt_extended import get_jwt_identity, jwt_required

import chemquest_website.elo as elo
from chemquest_website import app, db

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