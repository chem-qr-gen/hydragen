from flask import request
from flask_login import current_user, login_required

from chemquest_website import app, engine, meta

@app.route('/edit_profile', methods=['POST'])
@login_required
def edit_profile():
    # edits profile of user
    username = current_user.get_id()
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
