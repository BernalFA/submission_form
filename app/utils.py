import base64
import io
from io import BytesIO
from itertools import product

import openpyxl
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from app.config import ALLOWED_FIELDS


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
            if len(str(num)) == 1:
                plate_pos.append(f"{letter}0{num}")
            else:
                plate_pos.append(f"{letter}{num}")
        return plate_pos


def get_allowed_sql_fields():
    internal_fields = set(ALLOWED_FIELDS["internal"].values())
    external_fields = set(ALLOWED_FIELDS["external"].values())
    return internal_fields & external_fields


def make_input_valid(entry):
    numeric = ["p_num", "mw", "amount", "vol", "conc"]
    res = {}
    for key, value in entry.items():
        if key in get_allowed_sql_fields():
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


ALLOWED_EXTENSIONS = {"xlsx"}


def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


# created by asking ChatGPT for template validation
def validate_excel_template(file_bytes, collaborator):
    try:
        wb = openpyxl.load_workbook(BytesIO(file_bytes))
        sheet = wb.active  # or use wb["sheet1"] if specific sheet
        headers = [cell.value for cell in next(sheet.iter_rows(min_row=1, max_row=1))]

        return headers == list(ALLOWED_FIELDS[collaborator].keys())

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


def export_to_excel(user, compounds):
    order = list(ALLOWED_FIELDS[user.membership].values())
    col_names = list(ALLOWED_FIELDS[user.membership].keys())

    cmpd_data = [c.__dict__ for c in compounds]
    for row in cmpd_data:
        row.pop("_sa_instance_state", None)
        row.pop("id", None)
    df = pd.DataFrame(cmpd_data)
    df = df[order]
    df.columns = col_names
    filepath = f"/home/freddy/Documents/{user.username}_{user.membership}.xlsx"
    df.to_excel(filepath, index=False)


def rename_columns(df, membership):
    return df.rename(columns=ALLOWED_FIELDS[membership])
