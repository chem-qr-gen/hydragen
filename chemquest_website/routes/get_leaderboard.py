from flask import jsonify
from sqlalchemy import select
from chemquest_website import app, engine, meta
import json

@app.route('/get_leaderboard', methods=['GET'])
def get_leaderboard():
    # Return all users and their ratings
    users_table = meta.tables["users"]
    with engine.connect() as conn:
        users_by_rating = conn.execute(select(users_table.c.username, users_table.c.elo)).fetchall()
    users_by_rating = [user._asdict() for user in users_by_rating]

    return jsonify(users_by_rating)
