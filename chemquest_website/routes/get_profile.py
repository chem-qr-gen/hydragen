from flask_jwt_extended import get_jwt_identity, jwt_required
from chemquest_website import app, engine, meta

@app.route('/get_profile', methods=['GET'])
@jwt_required(optional=True)
def get_profile():
    # Returns all the profile information of the user
    username = get_jwt_identity()

    users_table = meta.tables["users"]
    with engine.connect() as conn:
        user_dict = conn.execute(
            users_table.select().where(users_table.c.username == username)
        ).fetchone()
    user_dict = user_dict._asdict()

    elo = user_dict["elo"] if "elo" in user_dict else 1000.0
    email = user_dict["email"]
    gender = user_dict["gender"]
    country = user_dict["country"]
    region = user_dict["region"]
    return {"username": username, "elo": elo, "email": email, "gender": gender, "country": country, "region": region}
