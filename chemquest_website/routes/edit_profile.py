from flask import request
from flask_jwt_extended import get_jwt_identity, jwt_required

from chemquest_website import app, engine, meta

@app.route('/edit_profile', methods=['POST'])
@jwt_required(optional=True)
def edit_profile():
    # edits profile of user
    username = get_jwt_identity()
    gender = request.json["gender"]
    country = request.json["country"]
    region = request.json["region"]

    users_table = meta.tables["users"]
    with engine.connect() as conn:
        conn.execute(users_table.update()
                     .where(users_table.c.username == username)
                     .values({"gender": gender, "country": country, "region": region}))
        conn.commit()

    return {"msg": "Profile updated successfully"}
