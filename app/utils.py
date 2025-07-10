import base64
import io
from io import BytesIO
from itertools import product

import openpyxl
from flask_wtf import FlaskForm
from rdkit import Chem
from rdkit.Chem import Draw
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


def smiles_to_png_base64(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except TypeError:
        return None
    img = Draw.MolToImage(mol, size=(200, 200))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()
