import pandas as pd
from flask import (
    Blueprint,
    Response,
    flash,
    redirect,
    render_template,
    request,
    send_from_directory,
    session,
    url_for,
)

from app.extensions import db
from app.forms import UserDataForm
from app.models import UserManager, CompoundManager
from app.utils import (
    PositionGenerator,
    allowed_file,
    # export_to_excel,
    make_input_valid,
    rename_columns,
    smiles_to_png_base64,
    validate_excel_template,
)


# define a generator for the plate position
position_generator = PositionGenerator()


# define a main blueprint
main_bp = Blueprint("main", __name__)


@main_bp.route("/", methods=["GET", "POST"])
def index():
    user_form = UserDataForm()
    if user_form.validate_on_submit():
        session["delivery"] = user_form.delivery.data
        session["sample_type"] = user_form.sample_type.data
        session["include_structures"] = user_form.include_structures.data

        affiliation = user_form.affiliation.data
        entry = UserManager(
            session_id=session["session_id"],
            username=user_form.username.data,
            email=user_form.email.data,
            affiliation=affiliation,
            delivery=user_form.delivery.data,
            include_structures=user_form.include_structures.data == "true"
        )
        try:
            db.session.add(entry)
            db.session.commit()
            return redirect(url_for("main.upload", affiliation=affiliation))
        except Exception as e:
            print(f"ERROR: {e}")
            return f"ERROR: {e}"
    return render_template("index.html", user_form=user_form)


@main_bp.route("/upload/<affiliation>", methods=["GET", "POST"])
def upload(affiliation):
    delivery = session["delivery"]
    sample_type = session["sample_type"]
    include_structures = session["include_structures"]
    compounds = CompoundManager.query.filter_by(session_id=session["session_id"]).all()
    if request.method == "POST":
        entry = request.form
        entry = make_input_valid(entry)
        if entry:
            entry["session_id"] = session["session_id"]
            if affiliation == "internal":
                entry["position"] = position_generator.get_position()
            elif affiliation == "external":
                entry["png"] = smiles_to_png_base64(entry["smiles"])
            new_entry = CompoundManager(**entry)
            try:
                db.session.add(new_entry)
                db.session.commit()
                return redirect(url_for("main.upload", affiliation=affiliation))
            except Exception as e:
                print(f"ERROR: {e}")
                return f"ERROR: {e}"
    if delivery == "vials_solid" and affiliation == "external":
        template = "upload_external_dry.html"
    elif affiliation == "internal":
        template = "upload_internal.html"
    else:
        template = "upload_external_solution.html"
    return render_template(template, compounds=compounds, delivery=delivery, sample_type=sample_type, include_structures=include_structures)


@main_bp.route("/summary/<affiliation>")
def summary(affiliation):
    compounds = CompoundManager.query.filter_by(session_id=session["session_id"]).all()
    user = UserManager.query.filter_by(session_id=session["session_id"]).all()[0]
    # export_to_excel(user, compounds)
    return render_template(f"summary_{affiliation}.html", user=user, compounds=compounds)


@main_bp.route("/end")
def end():
    return render_template("end.html")


@main_bp.route("/download/<affiliation>")
def download_template(affiliation):
    filename = f"{affiliation.upper()}_Compound_submission.xlsx"
    return send_from_directory("downloads", filename, as_attachment=True)


@main_bp.route("/upload_from_file/<affiliation>", methods=["GET", "POST"])
def upload_from_file(affiliation):
    compounds = CompoundManager.query.filter_by(session_id=session["session_id"]).all()
    if request.method == "POST":
        file = request.files["file"]

        if not file or file.filename == "":
            flash("No selected file")
            return redirect(request.url)

        if not allowed_file(file.filename):
            flash("Invalid file type!")
            return redirect(request.url)

        file_bytes = file.read()
        file.seek(0)  # Reset stream pointer after reading

        if not validate_excel_template(file_bytes, affiliation):
            flash("Uploaded Excel file does not match the expected template!")
            return redirect(request.url)

        df = pd.read_excel(file)
        df.dropna(
            subset=[col for col in df.columns if col != "Position"],
            how="all",
            inplace=True
        )
        df = rename_columns(df, affiliation)
        df["session_id"] = session["session_id"]
        if affiliation == "external":
            df["png"] = df["smiles"].apply(smiles_to_png_base64)

        for _, row in df.iterrows():
            new_entry = CompoundManager(**row)
            try:
                db.session.add(new_entry)
                db.session.commit()
            except Exception as e:
                return f"ERROR: {e}"

        return redirect(url_for("main.upload_from_file", affiliation=affiliation))
    return render_template(f"upload_{affiliation}_from_file.html", compounds=compounds)


@main_bp.route("/mol_png/<int:mol_id>")
def mol_png(mol_id):
    cmpd = CompoundManager.query.get_or_404(mol_id)
    if not cmpd.png:
        return Response("Image not found", status=404)
    html = f'<img src="data:image/png;base64,{cmpd.png}" width="200" height="200">'
    return html


@main_bp.route("/reset")
def reset_session():
    session_id = session.get("session_id")
    if session_id:
        UserManager.query.filter_by(session_id=session_id).delete()
        CompoundManager.query.filter_by(session_id=session_id).delete()
        db.session.commit()
    session.clear()
    return redirect(url_for("main.end"))
