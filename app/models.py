from app.extensions import db


# define model using SQLAlchemy
# based on YouTube video
# https://www.youtube.com/watch?v=45P3xQPaYxc
class CompoundManager(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    session_id = db.Column(db.String(100))
    position = db.Column(db.String, unique=True)

    # Shared
    stereo_comment = db.Column(db.String, nullable=True)
    mw = db.Column(db.String)
    amount = db.Column(db.Float)
    vol = db.Column(db.Float)
    conc = db.Column(db.Float)
    project = db.Column(db.String)
    comment = db.Column(db.String)

    # Internal-specific
    exp_name = db.Column(db.String, nullable=True)
    p_num = db.Column(db.String, nullable=True)

    # External-specific
    supplier = db.Column(db.String, nullable=True)
    supp_id = db.Column(db.String, nullable=True)
    producer = db.Column(db.String, nullable=True)
    trivial_name = db.Column(db.String, nullable=True)
    alt_name = db.Column(db.String, nullable=True)
    cas = db.Column(db.String, nullable=True)
    smiles = db.Column(db.String, nullable=True)
    annotation = db.Column(db.String, nullable=True)
    png = db.Column(db.Text, nullable=True)

    def __repr__(self):
        return f"Entry: {self.position}"


class UserManager(db.Model):
    __bind_key__ = "user"
    id = db.Column(db.Integer, primary_key=True)
    session_id = db.Column(db.String(100))
    username = db.Column(db.String)
    email = db.Column(db.String, unique=True)
    membership = db.Column(db.String)
    delivery = db.Column(db.String)
    include_structures = db.Column(db.Boolean)
