import os
from dataclasses import dataclass
from dotenv import load_dotenv
from typing import Union, Optional

# load key
load_dotenv()


# define export directory
EXPORT_DIR = "/home/freddy/Documents/"


# define Flask-SQLAlchemy configuration
class Config:
    SECRET_KEY = os.getenv("FLASK-WTFS_KEY")
    SQLALCHEMY_DATABASE_URI = "sqlite:///database.db"


# define column names configuration
@dataclass
class ColumnRule:
    required: bool
    db_name: str
    type: Optional[Union[str, int, float]] = None


# Define schemas
base_schema = {
    "Position": ColumnRule(required=True, db_name="position", type=str),
    "Stereo comment": ColumnRule(required=True, db_name="stereo_comment"),
    "Molecular weight": ColumnRule(required=False, db_name="mw"),
    "Amount (mg)": ColumnRule(required=True, db_name="amount"),
    "Volume (µl)": ColumnRule(required=True, db_name="vol"),
    "Conc. (mM)": ColumnRule(required=True, db_name="conc"),
    "Project name": ColumnRule(required=True, db_name="project", type=str),
    "Comment": ColumnRule(required=False, db_name="comment"),
}

internal_schema = {
    **base_schema,
    "Enso experiment name": ColumnRule(required=True, db_name="exp_name", type=str),
    "Product No": ColumnRule(required=False, db_name="p_num", type=int),
}

external_schema = {
    **base_schema,
    "Supplier": ColumnRule(required=False, db_name="supplier"),
    "Supplier ID": ColumnRule(required=True, db_name="supp_id"),
    "Producer": ColumnRule(required=False, db_name="producer"),
    "Trivial name": ColumnRule(required=False, db_name="trivial_name"),
    "Alternative names": ColumnRule(required=False, db_name="alt_name"),
    "CAS": ColumnRule(required=False, db_name="cas"),
    "SMILES": ColumnRule(required=False, db_name="smiles"),
    "Annotation": ColumnRule(required=False, db_name="annotation"),
}

ALLOWED_SCHEMA = {
    "internal": internal_schema,
    "external": external_schema
}
