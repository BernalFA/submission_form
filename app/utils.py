import base64
import io
from datetime import datetime
from itertools import product

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from app.config import ALLOWED_SCHEMA, ColumnRule, EXPORT_DIR


class PositionGenerator:
    """Utility for automatic assignment of a position in a 96-well plate"""

    def __init__(self):
        self._positions = iter(self._create_plate_positions())

    def get_position(self):
        """Assign next position in a 96-well plate if places still available from
        _create_plate_positions. If there is no place, no position will be assigned.

        Returns:
            str: plate position (e.g. B01)
        """
        try:
            return next(self._positions)
        except StopIteration:
            print("Plate full")

    def reset(self):
        """helper method to reset plate to first position (empty plate)"""
        self._positions = iter(self._create_plate_positions())

    def _create_plate_positions(self):
        """Generate basic 96-well plate naming for each well (position).

        Returns:
            list: positions in 96-well plate.
        """
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
    """Check the ALLOWED_SCHEMA and return the unique set of field names allowed in
    the database.

    Returns:
        set: allowed field names
    """
    internal_fields = [col.db_name for col in ALLOWED_SCHEMA["internal"].values()]
    external_fields = [col.db_name for col in ALLOWED_SCHEMA["external"].values()]
    return set(internal_fields + external_fields)


def make_input_valid(entry):
    """Manage user input and return it using the correct/expected type. It is necessary
    before commiting to database to avoid format errors.

    Args:
        entry (dict): user entry (typically Flask form/dict)

    Returns:
        dict: entries using the correct type.
    """
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
    """Simple check of file correctness by file extension.

    Args:
        filename (str): filename

    Returns:
        bool: True if file has correct extension.
    """
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def validate_excel_template(df, affiliation, delivery, sample_type, structures):
    """Helper function to validate Excel file submitted by user. It first customize
    the default schema based on user input. Then verification of expected columns and
    their content (when required) is carried out. Options must match those given by
    UserDataForm in forms.

    Args:
        df (pd.DataFrame): dataframe from user input Excel.
        affiliation (str): collaborator type.
        delivery (str): mode of sample delivery.
        sample_type (str): type of sample.
        structures (bool): whether chemical structures are shared.

    Returns:
        list: list of errors found (if any)
    """
    schema = customize_schema(affiliation, delivery, sample_type, structures)

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


def customize_schema(affiliation, delivery, sample_type, structures):
    """Change default schema based on user input (necessary to avoid conflicts upon
    file content validation). Options must match those given by UserDataForm in forms.

    Args:
        affiliation (str): collaborator type.
        delivery (str): mode of sample delivery.
        sample_type (str): type of sample.
        structures (bool): whether chemical structures are shared.

    Returns:
        dict: customized schema
    """
    schema = ALLOWED_SCHEMA[affiliation]
    if delivery == "vials_solid":
        schema["Volume (µl)"] = ColumnRule(required=False, db_name="vol")
        schema["Conc. (mM)"] = ColumnRule(required=False, db_name="conc")
    else:
        schema["Amount (mg)"] = ColumnRule(required=False, db_name="amount")

    if sample_type == "isolated":
        if affiliation == "external":
            if structures:
                schema["SMILES"] = ColumnRule(required=True, db_name="smiles")
            else:
                schema["Molecular weight"] = ColumnRule(required=True, db_name="mw")
    return schema


def smiles_to_png_base64(smiles):
    """Transform SMILES string into PNG image (compatible with HTML and Flask display).

    Args:
        smiles (str): SMILES string of specified compound.

    Returns:
        base64: enconding of PNG generated image for the molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except TypeError:
        return None
    img = Draw.MolToImage(mol, size=(200, 200))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()


def export_to_excel(user, compounds):
    """Create an Excel file containing the compounds the user has registered.

    Args:
        user (SQLAlchemy.Model.query): query of the UserManager table.
        compounds (SQLAlchemy.Model.query): query of the CompoundManager table.
    """
    order = [col.db_name for col in ALLOWED_SCHEMA[user.affiliation].values()]
    col_names = [col for col in ALLOWED_SCHEMA[user.affiliation].keys()]

    cmpd_data = [c.__dict__ for c in compounds]
    for row in cmpd_data:
        row.pop("_sa_instance_state", None)
        row.pop("id", None)
    df = pd.DataFrame(cmpd_data)
    df = df[order]
    df.columns = col_names

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filepath = f"{EXPORT_DIR}{user.username}_{user.affiliation}_{timestamp}.xlsx"

    with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
        df.to_excel(writer, index=False)


def rename_columns(df, affiliation):
    """Change names of columns in given dataframe to match those in the SQLAlchemy
    database.

    Args:
        df (pd.DataFrame): dataframe with user provided compound information.
        affiliation (str): type of collaboration.

    Returns:
        pd.DataFrame: dataframe with columns matching the SQL database fields.
    """
    rename_dict = {key: val.db_name for key, val in ALLOWED_SCHEMA[affiliation].items()}
    return df.rename(columns=rename_dict)
