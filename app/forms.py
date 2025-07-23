from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, RadioField, EmailField
from wtforms.validators import DataRequired


class UserDataForm(FlaskForm):
    username = StringField("Full name", validators=[DataRequired()])
    email = EmailField("Email", validators=[DataRequired()])
    membership = SelectField(
        "Affiliation", choices=[("internal", "MPI Dortmund"), ("external", "External")]
    )
    delivery = SelectField(
        "Samples delivered in", choices=[("vials", "Vials"), ("plate", "Plate")]
    )
    include_structures = RadioField(
        "Will you share the compounds' structures?",
        choices=[("true", "Yes"), ("false", "No")],
        default="true"
    )
