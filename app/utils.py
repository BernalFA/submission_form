from itertools import product

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
]


def create_plate_positions():
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


class PositionGenerator:
    def __init__(self):
        self._positions = iter(create_plate_positions())

    def get_position(self):
        try:
            return next(self._positions)
        except StopIteration:
            print("Plate full")


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
