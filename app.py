from flask import Flask, render_template, request, url_for
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chemquest.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Question(db.Model):
    id = db.Column(db.Integer, primary_key = True)
    text1 = db.Column(db.String(256))
    imgsrc = db.Column(db.String(128))
    text2 = db.Column(db.String(256))
    answers = db.Column(db.String(128), nullable = False)
    def __repr__(self):
        return f"<Question {self.id}>"
    def __iter__(self):
        for key in ("id", "text1", "imgsrc", "text2", "answers"):
            yield (key, getattr(self, key))

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/ms_questions')
def ms_questions():
    question_id = int(request.args.get('id'))
    return dict(Question.query.filter_by(id = question_id).first())

if __name__ == "__main__":
    app.run(debug = True, host = "0.0.0.0")