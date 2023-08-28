import tomllib

from flask import Flask
from flask_jwt_extended import JWTManager
from flask_mail import Mail
from flask_seasurf import SeaSurf
from sqlalchemy import create_engine, MetaData


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

# initialise the database engine
engine = create_engine(cfg["postgres_address"], isolation_level="AUTOCOMMIT")
meta = MetaData()
meta.reflect(bind=engine)

from chemquest_website.routes import *
