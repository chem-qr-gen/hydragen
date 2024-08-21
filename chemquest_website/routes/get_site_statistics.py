from sqlalchemy.orm import sessionmaker

from chemquest_website import app, engine, meta
from sqlalchemy import func

@app.route('/get_site_statistics')
def get_site_statistics():
    users_table = meta.tables["users"]
    attempts_table = meta.tables["attempts"]

    session = sessionmaker(bind=engine)
    session = session()
    users_count = session.query(func.count(users_table.c.username)).scalar()
    attempts_count = session.query(func.count(attempts_table.c.attempt_id)).scalar()

    return {"users_count": users_count, "attempts_count": attempts_count}
