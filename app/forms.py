from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, BooleanField, EmailField
from wtforms.validators import DataRequired, Optional


class UserDataForm(FlaskForm):
    username = StringField("Full name", validators=[DataRequired()])
    email = EmailField("Email", validators=[Optional()])
    affiliation = SelectField(
        "Affiliation", choices=[("internal", "MPI Dortmund"), ("external", "External")]
    )
    sample_type = SelectField(
        "Sample type",
        choices=[
            ("isolated", "Isolated compound"),
            ("mixture", "Mixture (e.g. extracts, fractions)")
        ],
    )
    delivery = SelectField(
        "Samples delivered in",
        choices=[
            ("vials_solid", "Vials (no solvent)"),
            ("vials_solution", "Vials (DMSO solution)"),
            ("plate", "Plate"),
        ]
    )
    include_structures = BooleanField(
        "Will you share the compounds' structures?",
        default="checked"
    )
