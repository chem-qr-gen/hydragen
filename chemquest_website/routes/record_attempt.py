import datetime
from flask import request
from flask_jwt_extended import get_jwt_identity, jwt_required

import chemquest_website.elo as elo
from chemquest_website import app, db

@app.route('/record_attempt', methods = ['POST'])
@jwt_required(optional = True)
def record_attempt():
    '''Record a question attempt made by the user, and relevant stats.'''
    hash_id = request.json['hashId']
    qid = request.json['qid']
    is_correct = request.json['isCorrect']

    jwt_identity = get_jwt_identity()
    if jwt_identity:

        # Check list of attempts to see if the question has already been attempted
        # If attempted, updates the attempt based on the number of wrong answers

        user_dict = db.users.find_one({"username": jwt_identity})
        question_dict = db.ms_data.find_one({"qid": qid})

        player_old_elo = user_dict["elo"] if "elo" in user_dict else 1000.0
        attempted_before = False

        if "attempts" not in user_dict:
            user_dict["attempts"] = []

        for i in range(len(user_dict["attempts"])):
            if user_dict["attempts"][i]["hash_id"] == hash_id:
                # This question has been attempted before. Modify the elo calculation and update player elo
                attempted_before = True                

                new_attempts = user_dict["attempts"]
                new_attempts[i]["is_correct"].append(is_correct)
                new_elo = elo.calculate_new_elo(
                    old_elo=new_attempts[i]["player_old_elo"],
                    question_elo=question_dict["difficulty"],
                    attempts=new_attempts[i]["is_correct"]
                )
                new_attempts[i]["player_new_elo"] = new_elo

                break
        
        if not attempted_before:
            # This question hasn't been attempted before. Add it to the attempts list and update the user
            new_elo = elo.calculate_new_elo(
                old_elo=player_old_elo,
                question_elo=question_dict["difficulty"],
                attempts=[is_correct]
            )
            this_attempt = {
                "timestamp": datetime.datetime.utcnow(),
                "hash_id": hash_id,
                "qid": qid,
                "is_correct": [is_correct],
                "player_old_elo": player_old_elo,
                "player_new_elo": new_elo,
            }
        
            new_attempts = user_dict["attempts"] if "attempts" in user_dict else []
            new_attempts.append(this_attempt)


        db.users.update_one(
            {"_id": user_dict["_id"]},
            {"$set": {
                "attempts": new_attempts,
                "elo": new_elo
            }},
            upsert = False)
    
        return {"msg": "Response recorded.", "new_elo": new_elo}
    return {"msg": "Not logged in, response not recorded."}