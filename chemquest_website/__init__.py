import tomllib

from flask import Flask
from flask_jwt_extended import JWTManager
from flask_mail import Mail
from flask_seasurf import SeaSurf
from FPSim2 import FPSim2Engine
from pymongo import MongoClient


app = Flask(__name__)

# load app configuration
with open("chemquest_website/chemquest_config.toml", "rb") as f:
    cfg = tomllib.load(f)
app.config.from_mapping(cfg)
#app.config['JWT_ACCESS_TOKEN_EXPIRES'] = datetime.timedelta(days = 1) # TODO: set up a JWT refresh token so the user is auto-logged out after inactivity.

# initialise the various addons and other services needed
jwt = JWTManager(app)
mail = Mail(app)
csrf = SeaSurf(app)
client = MongoClient(cfg["mongo_address"])
db = client.chemquest_db
fpe = FPSim2Engine("db/fingerprints.h5", in_memory_fps=False) # engine for molecule similarity search

from chemquest_website.routes import *
