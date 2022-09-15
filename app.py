import random
from flask import Flask, redirect, render_template, request, url_for
from pymongo import MongoClient

from mongo_address import mongo_address

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
client = MongoClient(mongo_address)
db = client.chemquest_db
ms_qns_collection = db.ms_qns

@app.route('/')
def index():
    '''Returns a skeleton index page for Flask. For the actual HTML, refer to files in static/src/views.'''
    return render_template("index.html")

@app.route('/ms_questions')
def ms_questions():
    '''Returns a question with the requested id, or a random id, in JSON form.'''
    question_id = request.args.get('id')
    if question_id == "random":
        return redirect(url_for('ms_questions', id = random.randint(1, 9)))
    question_response = ms_qns_collection.find_one({"qid": int(question_id)})
    question_response.pop("_id")
    return question_response

if __name__ == "__main__":
    app.run(debug = True, host = "0.0.0.0")
