import os
from dataclasses import dataclass
from dotenv import load_dotenv
from typing import Union, Optional

# load key
load_dotenv()


# define Flask-SQLAlchemy configuration
class Config:
    SECRET_KEY = os.getenv("FLASK-WTFS_KEY")
    SQLALCHEMY_DATABASE_URI = "sqlite:///database.db"


@dataclass
class ColumnRule:
    required: bool
    type: Optional[Union[str, int, float]] = None


base_schema = {
    "Position": ColumnRule(required=True, type=str),
    "Stereo comment": ColumnRule(required=False),
    "Molecular weight": ColumnRule(required=False),
    "Amount (mg)": ColumnRule(required=True),
    "Volume (µl)": ColumnRule(required=True),
    "Conc. (mM)": ColumnRule(required=True),
    "Project name": ColumnRule(required=True, type=str),
    "Comment": ColumnRule(required=False),
}

internal_schema = {
    **base_schema,
    "Enso experiment name": ColumnRule(required=True, type=str),
    "Product No": ColumnRule(required=False, type=int),
}

external_schema = {
    **base_schema,
    "Supplier": ColumnRule(required=False),
    "Supplier ID": ColumnRule(required=True),
    "Producer": ColumnRule(required=False),
    "Trivial name": ColumnRule(required=False),
    "Alternative names": ColumnRule(required=False),
    "CAS": ColumnRule(required=False),
    "SMILES": ColumnRule(required=False),
    "Annotation": ColumnRule(required=False),
}

ALLOWED_SCHEMA = {
    "internal": internal_schema,
    "external": external_schema
}


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
