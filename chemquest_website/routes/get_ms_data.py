import random
from flask import redirect, request, url_for

from chemquest_website import app, db

@app.route('/ms_questions_new')
def ms_questions_new():
    '''Returns a set of MS data based on the molecule id provided, or a random one if id == "random".'''
    question_id = request.args.get('id')
    if question_id == "random":
        return redirect(url_for('ms_questions_new', id = random.randint(0, db.ms_data.count_documents({}) - 1)))
    question_response = db.ms_data.find_one({"qid": int(question_id)})
    question_response.pop("_id") # remove the internal mongodb id from the returned JSON
    return question_response