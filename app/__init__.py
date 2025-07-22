import uuid
from flask import Flask, session

from app.config import Config
from app.extensions import db
from app.views import main_bp


# define application factory
def create_app():
    # initialize app and configure it
    app = Flask(__name__)
    app.config.from_object(Config)
    # initialize extensions (SQLAlchemy)
    db.init_app(app)

    # Define session logic (session id)
    @app.before_request
    def assign_session():
        if "session_id" not in session:
            session["session_id"] = str(uuid.uuid4())

    # add blueprint
    app.register_blueprint(main_bp)

    return app
