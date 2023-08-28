import random

from flask import redirect, request, url_for
from sqlalchemy import func, select

from chemquest_website import app, engine, meta

@app.route('/ms_questions_new')
def ms_questions_new():
    '''Returns a set of MS data based on the molecule id provided, or a random one if id == "random".'''
    question_id = request.args.get('id')

    ms_data_table = meta.tables["ms_data"]

    # pick a random question from the table
    if question_id == "random":
        with engine.connect() as conn:
            max_qid = conn.execute(select(func.max(ms_data_table.c.qid))).fetchone()[0]
        return redirect(url_for('ms_questions_new', id = random.randint(0, max_qid)))
    
    ms_data_table = meta.tables["ms_data"]
    with engine.connect() as conn:
        question_row = conn.execute(
            ms_data_table.select().where(ms_data_table.c.qid == int(question_id))
        ).fetchone()
        question_response = dict(question_row._mapping)

    # question_response = db.ms_data.find_one({"qid": int(question_id)})
    # question_response.pop("_id") # remove the internal mongodb id from the returned JSON
    return question_response