from app.extensions import db


# define model using SQLAlchemy
# based on YouTube video
# https://www.youtube.com/watch?v=45P3xQPaYxc
class CompoundManagerInternal(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    position = db.Column(db.String, unique=True)
    exp_name = db.Column(db.String)
    stereo_comment = db.Column(db.String)
    # p_num = db.Column(db.Integer)
    # mw = db.Column(db.Float)
    p_num = db.Column(db.String)
    mw = db.Column(db.String)
    amount = db.Column(db.Float)
    vol = db.Column(db.Float)
    conc = db.Column(db.Float)
    project = db.Column(db.String)
    comment = db.Column(db.String)

    def __repr__(self):
        return f"Entry: {self.position}"


class CompoundManagerExternal(db.Model):
    id = db.Column(db.Integer, primary_key=True)


class UserManager(db.Model):
    __bind_key__ = "user"
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String)
    membership = db.Column(db.String)
