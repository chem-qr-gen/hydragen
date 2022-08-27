import random
from flask import Flask, redirect, render_template, request, url_for
from flask_sqlalchemy import SQLAlchemy
from init_vars import init_vars

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chemquest.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Question(db.Model):
    '''A mass spectrometry question. For use with SQLAlchemy.'''
    id = db.Column(db.Integer, primary_key = True)
    text1 = db.Column(db.String(256))
    imgsrc = db.Column(db.String(128))
    text2 = db.Column(db.String(256))
    hint1 = db.Column(db.String(256))
    hint2 = db.Column(db.String(256))
    answers = db.Column(db.String(128), nullable = False)
    
    @init_vars
    def __init__(self, id, text1, imgsrc, text2, hint1, hint2, answers):
        pass
    
    def __repr__(self):
        return f"<Question {self.id}>"
    
    def __iter__(self):
        for key in ("id", "text1", "imgsrc", "text2", "hint1", "hint2", "answers"):
            yield (key, getattr(self, key))

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
    return dict(Question.query.filter_by(id = int(question_id)).first())

if __name__ == "__main__":
    app.run(debug = True, host = "0.0.0.0")
