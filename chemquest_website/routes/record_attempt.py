import datetime
from flask import request
from flask_login import current_user
from sqlalchemy.dialects.postgresql import insert

import chemquest_website.elo as elo
from chemquest_website import app, engine, meta

# TODO: add number of hints to the calculation
@app.route('/record_attempt', methods = ['POST'])
def record_attempt():
    '''Record a question attempt made by the user, and relevant stats.'''
    hash_id = request.json['hashId']
    qid = request.json['qid']
    is_correct = request.json['isCorrect']
    generateTimestamp = request.json["generate_timestamp"]
    user_ip = request.remote_addr

    users_table = meta.tables["users"]
    ms_data_table = meta.tables["ms_data"]
    attempts_table = meta.tables["attempts"]

    if current_user.is_authenticated:
        
        # Check list of attempts to see if the question has already been attempted
        # If attempted, updates the attempt based on the number of wrong answers

        # Retrieve the user and question from the database
        with engine.connect() as conn:
            user_dict = conn.execute(
                users_table.select().where(users_table.c.username == current_user.get_id())
            ).fetchone()
            question_dict = conn.execute(
                ms_data_table.select().where(ms_data_table.c.qid == qid)
            ).fetchone()

        if user_dict:
            user_dict = user_dict._mapping

        if question_dict:
            question_dict = question_dict._mapping

        player_old_elo = user_dict["elo"] if "elo" in user_dict else 1000.0

        # check if the current question has been attempted before
        with engine.connect() as conn:
            attempt = conn.execute(
                attempts_table.select().where(attempts_table.c.attempt_id == hash_id)
            ).fetchone()

        if attempt:
            attempt = dict(attempt._mapping)
            attempt["is_correct"].append(is_correct)
            new_elo = elo.calculate_new_elo(
                old_elo=attempt["player_old_elo"],
                question_elo=question_dict["difficulty"],
                attempts=attempt["is_correct"]
            )
            attempt["player_new_elo"] = new_elo

        else:
            # This question hasn't been attempted before. Add it to the attempts list and update the user
            new_elo = elo.calculate_new_elo(
                old_elo=player_old_elo,
                question_elo=question_dict["difficulty"],
                attempts=[is_correct]
            )
            attempt = {
                "attempt_id": hash_id,
                "username": current_user.get_id(),
                "qid": qid,
                "generate_timestamp": generateTimestamp,
                "timestamp": datetime.datetime.utcnow(),
                "is_correct": [is_correct],
                "player_old_elo": player_old_elo,
                "player_new_elo": new_elo,
                "ip_address": user_ip,
            }

        with engine.connect() as conn:
            conn.execute(
                insert(attempts_table).values(**attempt).on_conflict_do_update(
                    index_elements=["attempt_id"],
                    set_={
                        "is_correct": attempt["is_correct"],
                        "player_new_elo": attempt["player_new_elo"],
                    }
                )
            )
    
        return {"msg": "Response recorded.", "new_elo": new_elo}
    
    # User not logged in, record attempt without username
    else:
        # Retrieve the question from the database
        with engine.connect() as conn:
            question_dict = conn.execute(
                ms_data_table.select().where(ms_data_table.c.qid == qid)
            ).fetchone()

        if question_dict:
            question_dict = question_dict._mapping

        # check if the current question has been attempted before
        with engine.connect() as conn:
            attempt = conn.execute(
                attempts_table.select().where(attempts_table.c.attempt_id == hash_id)
            ).fetchone()
        
        if attempt:
            attempt = dict(attempt._mapping)
            attempt["is_correct"].append(is_correct)

        else:
            attempt = {
                "attempt_id": hash_id,
                "qid": qid,
                "generate_timestamp": generateTimestamp,
                "timestamp": datetime.datetime.utcnow(),
                "is_correct": [is_correct],
                "ip_address": user_ip,
            }

        with engine.connect() as conn:
            conn.execute(
                insert(attempts_table).values(**attempt).on_conflict_do_update(
                    index_elements=["attempt_id"],
                    set_={
                        "is_correct": attempt["is_correct"],
                    }
                )
            )

    return {"msg": "Response recorded. Not logged in."}
