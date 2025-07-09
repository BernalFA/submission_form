from io import BytesIO
from itertools import product

import openpyxl
# import pandas as pd
from flask_wtf import FlaskForm
from wtforms import StringField, SelectField
from wtforms.validators import DataRequired


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


EXPECTED_HEADERS = {
    "internal": [
        "Position",
        "Enso experiment name",
        "Stereo comment",
        "Product No",
        "Molecular weight",
        "Amount (mg)",
        "Volume (µl)",
        "Conc. (mM)",
        "Project name",
        "Comment",
    ],
    "external": [
        "Position",
        "Supplier",
        "Supplier ID",
        "Producer",
        "Stereo comment",
        "Molecular weight",
        "Amount (mg)",
        "Volume (µl)",
        "Conc. (mM)",
        "Project name",
        "Trivial name",
        "Alternative names",
        "CAS",
        "SMILES",
        "Annotation",
        "Comment",
    ],
}


# created by asking ChatGPT for template validation
def validate_excel_template(file_bytes, collaborator):
    try:
        wb = openpyxl.load_workbook(BytesIO(file_bytes))
        sheet = wb.active  # or use wb["sheet1"] if specific sheet
        headers = [cell.value for cell in next(sheet.iter_rows(min_row=1, max_row=1))]

        return headers == EXPECTED_HEADERS[collaborator]

    except Exception as e:
        print("Excel validation error:", e)
        return False

# SCHEMA = {
#     "internal": {
#         "Position": {"required": True, "type": str},
#         "Enso experiment name": {"required": True, "type": str},
#         "Stereo comment": {"required": False, "type": str},
#         "Product No": {"required": False, "type": int},
#         "Molecular weight": {"required": False, "type": float},
#         "Amount (mg)": {"required": True, "type": float},
#         "Volume (µl)": {"required": True, "type": float},
#         "Conc. (mM)": {"required": True, "type": float},
#         "Project name": {"required": True, "type": str},
#         "Comment": {"required": False, "type": str},
#     },
#     "external": {
#         "Position": {"required": True, "type": str},
#         "Supplier": {"required": False, "type": str},
#         "Supplier ID": {"required": False, "type": str},
#         "Producer": {"required": False, "type": str},
#         "Stereo comment": {"required": False, "type": str},
#         "Molecular weight": {"required": False, "type": float},
#         "Amount (mg)": {"required": True, "type": float},
#         "Volume (µl)": {"required": True, "type": float},
#         "Conc. (mM)": {"required": True, "type": float},
#         "Project name": {"required": True, "type": str},
#         "Trivial name": {"required": False, "type": str},
#         "Alternative names": {"required": False, "type": str},
#         "CAS": {"required": False, "type": str},
#         "SMILES": {"required": True, "type": str},
#         "Annotation": {"required": False, "type": str},
#         "Comment": {"required": False, "type": str},
#     },
# }

# def validate_excel_template(df, collaborator):
#     errors = []

#     # check for missing or extra columns
#     expected_columns = list(SCHEMA[collaborator].keys())
#     if set(df.columns) != set(expected_columns):
#         missing = set(expected_columns) - set(df.columns)
#         extra = set(df.columns) - set(expected_columns)
#         if missing:
#             errors.append(f"Missing columns: {', '.join(missing)}")
#         if extra:
#             errors.append(f"Unexpected columns: {', '.join(extra)}")

#     # validate each column according to schema
#     for col, rules in SCHEMA.items():
#         if rules.get("required") and df[col].isnull().any():
#             errors.append(f"Missing data in required column: {col}")

#         expected_type = rules.get("type")
#         if expected_type:
#             if not df[col].map(
# lambda x: isinstance(x, expected_type) or pd.isnull(x)
# ).all():
#                 errors.append(
# f"Invalid type in column '{col}': expected {expected_type.__name__}"
# )

#     return errors
