from flask import Blueprint, render_template, request, redirect, url_for

from app.extensions import db
from app.models import UserManager, CompoundManagerInternal, CompoundManagerExternal
from app.utils import PositionGenerator, make_input_valid, UserDataForm


# define a generator for the plate position
position_generator = PositionGenerator()


# define a main blueprint
main_bp = Blueprint("main", __name__)


@main_bp.route("/", methods=["GET", "POST"])
def index():
    user_form = UserDataForm()
    if user_form.validate_on_submit():
        entry = UserManager(username=user_form.username.data,
                            membership=user_form.membership.data)
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


@main_bp.route("/upload_internal", methods=["GET", "POST"])
def upload_internal():
    compounds = CompoundManagerInternal.query.all()
    if request.method == "POST":
        entry = request.form
        entry = make_input_valid(entry)
        if entry:
            entry["position"] = position_generator.get_position()
            new_entry = CompoundManagerInternal(**entry)
            try:
                db.session.add(new_entry)
                db.session.commit()
                return redirect(url_for("main.upload_internal"))
            except Exception as e:
                print(f"ERROR: {e}")
                return f"ERROR: {e}"
    return render_template("upload_internal.html", compounds=compounds)


@main_bp.route("/upload_external")
def upload_external():
    return render_template("upload_external.html")


@main_bp.route("/summary_internal")
def summary_internal():
    compounds = CompoundManagerInternal.query.all()
    user = UserManager.query.all()[0]
    return render_template("summary_internal.html", user=user, compounds=compounds)


@main_bp.route("/summary_external")
def summary_external():
    compounds = CompoundManagerExternal.query.all()
    user = UserManager.query.all()[0]
    return render_template("summary_external.html", user=user, compounds=compounds)


@main_bp.route("/end")
def end():
    return render_template("end.html")
