from flask import Flask

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
    # add blueprint
    app.register_blueprint(main_bp)

    return app
