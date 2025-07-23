import os
from dotenv import load_dotenv

# load key
load_dotenv()


# define Flask-SQLAlchemy configuration
class Config:
    SECRET_KEY = os.getenv("FLASK-WTFS_KEY")
    SQLALCHEMY_DATABASE_URI = "sqlite:///compounds.db"
    SQLALCHEMY_BINDS = {
        "user": "sqlite:///user.db",
    }


# Define fields for database creation
# Excel name : html/SQL name
ALLOWED_FIELDS = {
    "internal": {
        "Position": "position",
        "Enso experiment name": "exp_name",
        "Stereo comment": "stereo_comment",
        "Product No": "p_num",
        "Molecular weight": "mw",
        "Amount (mg)": "amount",
        "Volume (µl)": "vol",
        "Conc. (mM)": "conc",
        "Project name": "project",
        "Comment": "comment",
    },
    "external": {
        "Position": "position",
        "Supplier": "supplier",
        "Supplier ID": "supp_id",
        "Producer": "producer",
        "Stereo comment": "stereo_comment",
        "Molecular weight": "mw",
        "Amount (mg)": "amount",
        "Volume (µl)": "vol",
        "Conc. (mM)": "conc",
        "Project name": "project",
        "Trivial name": "trivial_name",
        "Alternative names": "alt_name",
        "CAS": "cas",
        "SMILES": "smiles",
        "Annotation": "annotation",
        "Comment": "comment",
    },
}
