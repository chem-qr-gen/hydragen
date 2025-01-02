import tomllib

from flask import Flask
from flask_mail import Mail
from flask_login import LoginManager, UserMixin
from sqlalchemy import create_engine, MetaData

app = Flask(__name__)

# load app configuration
with open("chemquest_website/chemquest_config.toml", "rb") as f:
    cfg = tomllib.load(f)
app.config.from_mapping(cfg)

# initialise the various addons and other services needed
mail = Mail(app)
login_manager = LoginManager(app)

# initialise the database engine
engine = create_engine(cfg["postgres_address"], isolation_level="AUTOCOMMIT", pool_pre_ping=True)
meta = MetaData()
meta.reflect(bind=engine)


#flask-login user loader
class User(UserMixin):
    def __init__(self, username):
        self.username = username

    def get_id(self):
        return self.username

@login_manager.user_loader
def load_user(user_id):
    return User(user_id) # TODO: figure out what this actually needs to be restricted to


from chemquest_website.routes import *
