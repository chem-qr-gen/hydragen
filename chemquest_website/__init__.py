from flask import Flask
from flask_jwt_extended import JWTManager
from flask_mail import Mail
from flask_seasurf import SeaSurf
from FPSim2 import FPSim2Engine
from pymongo import MongoClient

from chemquest_website.chemquest_secrets import * # secret addresses, change or remove for local testing

app = Flask(__name__)

# app configuration
app.config['JSON_SORT_KEYS'] = False # don't sort keys in alphabetical order when sending JSON responses
app.config['SECRET_KEY'] = csrf_secret
app.config['JWT_SECRET_KEY'] = jwt_secret
app.config['MAIL_SERVER'] = mail_server_secret
app.config['MAIL_DEFAULT_SENDER'] = mail_sender_secret
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USERNAME'] = mail_username_secret
app.config['MAIL_PASSWORD'] = mail_password_secret
#app.config['JWT_ACCESS_TOKEN_EXPIRES'] = datetime.timedelta(days = 1) # TODO: set up a JWT refresh token so the user is auto-logged out after inactivity.

# initalise the various addons and other services needed
jwt = JWTManager(app)
mail = Mail(app)
csrf = SeaSurf(app)
client = MongoClient(mongo_address)
db = client.chemquest_db
fpe = FPSim2Engine("db/fingerprints.h5", in_memory_fps=False) # engine for molecule similarity search

from chemquest_website.routes import *
