from flask import Flask, render_template, request, redirect, url_for
from flask_sqlalchemy import SQLAlchemy
from flask_wtf import FlaskForm
from wtforms import StringField, SelectField  #, IntegerField, FloatField, SubmitField
from wtforms.validators import DataRequired

from utils import PositionGenerator, make_input_valid


class Config:
    SECRET_KEY = "you-will-never-guess"


# application set up
app = Flask(__name__)

app.config.from_object(Config)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///compounds.db"
app.config["SQLALCHEMY_BINDS"] = {
    "user": "sqlite:///user.db"
}
db = SQLAlchemy(app)


# define model using SQLAlchemy
# based on YouTube video
# https://www.youtube.com/watch?v=45P3xQPaYxc
class DataManager(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    position = db.Column(db.String, unique=True)
    exp_name = db.Column(db.String)
    stereo_comment = db.Column(db.String)
    # p_num = db.Column(db.Integer)
    # mw = db.Column(db.Float)
    p_num = db.Column(db.String)
    mw = db.Column(db.String)
    amount = db.Column(db.Float)
    vol = db.Column(db.Float)
    conc = db.Column(db.Float)
    project = db.Column(db.String)
    comment = db.Column(db.String)

    def __repr__(self):
        return f"Entry: {self.position}"


class UserManager(db.Model):
    __bind_key__ = "user"
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String)
    membership = db.Column(db.String)


class UserDataForm(FlaskForm):
    username = StringField("Full name", validators=[DataRequired()])
    membership = SelectField(
        "Affiliation", choices=[("mpi_do", "MPI Dortmund"), ("external", "External")]
    )
    # submit = SubmitField("Submit")


# class InternalCompoundForm(FlaskForm):
#     exp_name = StringField("Enso experiment name", validators=[DataRequired()])
#     stereo_comment = StringField("Stereo comment")
#     p_num = IntegerField("Product No.")
#     mw = FloatField("Molecular Weight")
#     amount = FloatField("Amount (mg)", validators=[DataRequired()])
#     vol = FloatField("Volume (uL)", validators=[DataRequired()])
#     conc = FloatField("Conc. (mM)", validators=[DataRequired()])
#     project = StringField("Project Name", validators=[DataRequired()])
#     comment = StringField("Comment")

#     @classmethod
#     def new(cls):
#         form = cls()
#         form.exp_name.data = ""
#         return form


position_generator = PositionGenerator()


@app.route("/", methods=["GET", "POST"])
def index():
    user_form = UserDataForm()
    if user_form.validate_on_submit():
        entry = UserManager(username=user_form.username.data,
                            membership=user_form.membership.data)
        print(0, user_form.username.data)
        try:
            db.session.add(entry)
            db.session.commit()
            if user_form.membership.data == "mpi_do":
                return redirect("/upload_internal")
            else:
                return redirect("/upload_external")
        except Exception as e:
            print(f"ERROR: {e}")
            return f"ERROR: {e}"
    return render_template("index.html", user_form=user_form)


@app.route("/upload_internal", methods=["GET", "POST"])
def upload_internal():
    compounds = DataManager.query.all()
    if request.method == "POST":
    # cmp_form = InternalCompoundForm.new()
    # if cmp_form.validate_on_submit():
    #     print(cmp_form.data)
        entry = request.form
        # print(entry)
        entry = make_input_valid(entry)
        if entry:
            entry["position"] = position_generator.get_position()
            new_entry = DataManager(**entry)
            try:
                db.session.add(new_entry)
                db.session.commit()
                return redirect(url_for("upload_internal"))
            except Exception as e:
                print(f"ERROR: {e}")
                return f"ERROR: {e}"
    return render_template("upload_internal.html", compounds=compounds)


@app.route("/upload_external")
def upload_external():
    return render_template("upload_external.html")


@app.route("/summary")
def summary():
    compounds = DataManager.query.all()
    user = UserManager.query.all()[0]
    return render_template("summary.html", user=user, compounds=compounds)


@app.route("/end")
def end():
    return render_template("end.html")


if __name__ == "__main__":
    with app.app_context():
        db.drop_all()
        db.create_all()

        app.run(debug=True)
