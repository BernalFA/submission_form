from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, RadioField, EmailField
from wtforms.validators import DataRequired, Email


class UserDataForm(FlaskForm):
    username = StringField("Full name", validators=[DataRequired()])
    email = EmailField("Email", validators=[
        DataRequired(), Email("Invalid email address")
    ])
    membership = SelectField(
        "Affiliation", choices=[("internal", "MPI Dortmund"), ("external", "External")]
    )
    delivery = SelectField(
        "Compounds delivered in", choices=[("vials", "Vials"), ("plate", "Plate")]
    )
    include_structures = RadioField(
        "Will you share the compounds' structures?",
        choices=[("true", "Yes"), ("false", "No")],
        default="true"
    )
