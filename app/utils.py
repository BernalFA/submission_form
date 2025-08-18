import base64
import io
from itertools import product

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from app.config import ALLOWED_SCHEMA


class PositionGenerator:
    def __init__(self):
        self._positions = iter(self._create_plate_positions())

    def get_position(self):
        try:
            return next(self._positions)
        except StopIteration:
            print("Plate full")

    def reset(self):
        self._positions = iter(self._create_plate_positions())

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
    internal_fields = [col.db_name for col in ALLOWED_SCHEMA["internal"].values()]
    external_fields = [col.db_name for col in ALLOWED_SCHEMA["external"].values()]
    return set(internal_fields + external_fields)


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
    order = [col.db_name for col in ALLOWED_SCHEMA[user.affiliation].values()]
    col_names = [col.db_name for col in ALLOWED_SCHEMA[user.affiliation].keys()]

    cmpd_data = [c.__dict__ for c in compounds]
    for row in cmpd_data:
        row.pop("_sa_instance_state", None)
        row.pop("id", None)
    df = pd.DataFrame(cmpd_data)
    df = df[order]
    df.columns = col_names
    filepath = f"/home/freddy/Documents/{user.username}_{user.affiliation}.xlsx"
    df.to_excel(filepath, index=False)


def rename_columns(df, affiliation):
    rename_dict = {key: val.db_name for key, val in ALLOWED_SCHEMA[affiliation].items()}
    return df.rename(columns=rename_dict)
