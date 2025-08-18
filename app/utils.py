from itertools import product

import pandas as pd
from flask_wtf import FlaskForm
from wtforms import StringField, SelectField
from wtforms.validators import DataRequired

from app.config import ALLOWED_SCHEMA


ALLOWED_FIELDS = [
    "exp_name",
    "stereo_comment",
    "p_num",
    "mw",
    "amount",
    "vol",
    "conc",
    "project",
    "comment",
    "position",
    "supplier",
    "supp_id",
    "producer",
    "trivial_name",
    "alt_name",
    "cas",
    "smiles",
    "annotation"
]


class PositionGenerator:
    def __init__(self):
        self._positions = iter(self._create_plate_positions())

    def get_position(self):
        try:
            return next(self._positions)
        except StopIteration:
            print("Plate full")

    def _create_plate_positions(self):
        nums = range(1, 13)
        letters = [chr(character) for character in range(ord("A"), ord("H") + 1)]

        plate_pos = []
        for num, letter in product(nums, letters):
            # print(letter, num)
            if len(str(num)) == 1:
                plate_pos.append(f"{letter}0{num}")
            else:
                plate_pos.append(f"{letter}{num}")
        return plate_pos


def make_input_valid(entry):
    numeric = ["p_num", "mw", "amount", "vol", "conc"]
    res = {}
    for key, value in entry.items():
        if key in ALLOWED_FIELDS:
            if key in numeric:
                if key == "p_num":
                    try:
                        res[key] = int(value)
                    except ValueError:
                        res[key] = value
                else:
                    try:
                        res[key] = float(value)
                    except ValueError:
                        res[key] = value
            else:
                res[key] = value
    return res


class UserDataForm(FlaskForm):
    username = StringField("Full name", validators=[DataRequired()])
    membership = SelectField(
        "Affiliation", choices=[("mpi_do", "MPI Dortmund"), ("external", "External")]
    )


ALLOWED_EXTENSIONS = {"xlsx"}


def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def validate_excel_template(df, collaborator):
    schema = ALLOWED_SCHEMA[collaborator]

    errors = []

    # check for missing or extra columns
    expected_columns = list(schema.keys())
    missing = set(expected_columns) - set(df.columns)
    extra = set(df.columns) - set(expected_columns)
    if missing:
        errors.append(f"Missing columns: {', '.join(missing)}")
    if extra:
        errors.append(f"Unexpected columns: {', '.join(extra)}")

    # validate each column according to schema
    for col, rules in schema.items():
        if rules.required:
            if df[col].isnull().any():
                errors.append(f"Missing data in required column: {col}")

            if rules.type is not None:
                if not df[col].map(
                    lambda x: isinstance(x, rules.type) or pd.isnull(x)
                ).all():
                    errors.append(
                        f"Invalid type in column '{col}': expected {rules.type}"
                    )

    return errors
