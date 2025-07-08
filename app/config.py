import os
from dotenv import load_dotenv

# load key
load_dotenv()


# define Flask-SQLAlchemy configuration
class Config:
    SECRET_KEY = os.getenv("FLASK-WTFS_KEY")
    SQLALCHEMY_DATABASE_URI = "sqlite:///compounds_internal.db"
    SQLALCHEMY_BINDS = {
        "user": "sqlite:///user.db",
        "compounds_external": "sqlite:///compounds_external.db"
    }
