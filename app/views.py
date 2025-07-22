import pandas as pd
from flask import Blueprint, render_template, request, redirect, url_for, send_from_directory, flash
# from werkzeug.utils import secure_filename

from app.extensions import db
from app.models import UserManager, CompoundManagerInternal, CompoundManagerExternal
from app.utils import PositionGenerator, make_input_valid, UserDataForm, allowed_file, validate_excel_template, export_to_excel


# define a generator for the plate position
position_generator = PositionGenerator()


# define a main blueprint
main_bp = Blueprint("main", __name__)


@main_bp.route("/", methods=["GET", "POST"])
def index():
    user_form = UserDataForm()
    if user_form.validate_on_submit():
        membership = user_form.membership.data
        entry = UserManager(username=user_form.username.data,
                            membership=membership)
        try:
            db.session.add(entry)
            db.session.commit()
            return redirect(url_for("main.upload", membership=membership))
        except Exception as e:
            print(f"ERROR: {e}")
            return f"ERROR: {e}"
    return render_template("index.html", user_form=user_form)


@main_bp.route("/upload/<membership>", methods=["GET", "POST"])
def upload(membership):
    if membership == "internal":
        model = CompoundManagerInternal
        template = "upload_internal.html"
    elif membership == "external":
        model = CompoundManagerExternal
        template = "upload_external.html"

    compounds = model.query.all()
    if request.method == "POST":
        entry = request.form
        entry = make_input_valid(entry)
        if entry:
            if membership == "internal":
                entry["position"] = position_generator.get_position()
            new_entry = model(**entry)
            try:
                db.session.add(new_entry)
                db.session.commit()
                return redirect(url_for("main.upload", membership=membership))
            except Exception as e:
                print(f"ERROR: {e}")
                return f"ERROR: {e}"
    return render_template(template, compounds=compounds)


@main_bp.route("/summary_internal")
def summary_internal():
    compounds = CompoundManagerInternal.query.all()
    user = UserManager.query.all()[0]
    export_to_excel(user, compounds)
    return render_template("summary_internal.html", user=user, compounds=compounds)


@main_bp.route("/summary_external")
def summary_external():
    compounds = CompoundManagerExternal.query.all()
    user = UserManager.query.all()[0]
    export_to_excel(user, compounds)
    return render_template("summary_external.html", user=user, compounds=compounds)


@main_bp.route("/end")
def end():
    return render_template("end.html")


@main_bp.route("/download_internal")
def download_internal():
    filename = "INTERNAL_Compound_submission.xlsx"
    return send_from_directory("downloads", filename, as_attachment=True)


@main_bp.route("/download_external")
def download_external():
    filename = "EXTERNAL_collaboration_Compound_submission.xlsx"
    return send_from_directory("downloads", filename, as_attachment=True)


@main_bp.route("/upload_internal_from_file", methods=["GET", "POST"])
def upload_internal_from_file():
    compounds = CompoundManagerInternal.query.all()
    if request.method == "POST":
        file = request.files["file"]

        # if "file" not in request.files:
        #     flash("No file part")
        #     return redirect(request.url)

        if not file or file.filename == "":
            flash("No selected file")
            return redirect(request.url)

        if not allowed_file(file.filename):
            flash("Invalid file type!")
            return redirect(request.url)

        # filename = secure_filename(file.filename)
        file_bytes = file.read()
        file.seek(0)  # Reset stream pointer after reading

        if not validate_excel_template(file_bytes, "internal"):
            flash("Uploaded Excel file does not match the expected template!")
            return redirect(request.url)

        # filepath = os.path.join(main_bp.root_path, "uploads", filename)
        # file.save(filepath)
        # flash("File uploaded successfully!")
        df = pd.read_excel(file)
        df.dropna(subset=[col for col in df.columns if col != "Position"], how="all", inplace=True)

        for _, row in df.iterrows():
            new_entry = CompoundManagerInternal(
                position=row["Position"],
                exp_name=row["Enso experiment name"],
                stereo_comment=row["Stereo comment"],
                p_num=row["Product No"],
                mw=row["Molecular weight"],
                amount=row["Amount (mg)"],
                vol=row["Volume (µl)"],
                conc=row["Conc. (mM)"],
                project=row["Project name"],
                comment=row["Comment"],
            )
            try:
                db.session.add(new_entry)
                db.session.commit()
            except Exception as e:
                return f"ERROR: {e}"

        return redirect(url_for("main.upload_internal_from_file"))
    return render_template("upload_internal_from_file.html", compounds=compounds)


@main_bp.route("/upload_external_from_file", methods=["GET", "POST"])
def upload_external_from_file():
    compounds = CompoundManagerExternal.query.all()
    if request.method == "POST":
        file = request.files["file"]

        if not file or file.filename == "":
            flash("No selected file")
            return redirect(request.url)

        if not allowed_file(file.filename):
            flash("Invalid file type!")
            return redirect(request.url)

        # filename = secure_filename(file.filename)
        file_bytes = file.read()
        file.seek(0)  # Reset stream pointer after reading

        if not validate_excel_template(file_bytes, "external"):
            flash("Uploaded Excel file does not match the expected template!")
            return redirect(request.url)

        # filepath = os.path.join(main_bp.root_path, "uploads", filename)
        # file.save(filepath)
        # flash("File uploaded successfully!")
        df = pd.read_excel(file)
        df.dropna(subset=[col for col in df.columns if col != "Position"], how="all", inplace=True)

        for _, row in df.iterrows():
            new_entry = CompoundManagerExternal(
                position=row["Position"],
                supplier=row["Supplier"],
                supp_id=row["Supplier ID"],
                producer=row["Producer"],
                stereo_comment=row["Stereo comment"],
                mw=row["Molecular weight"],
                amount=row["Amount (mg)"],
                vol=row["Volume (µl)"],
                conc=row["Conc. (mM)"],
                project=row["Project name"],
                trivial_name=row["Trivial name"],
                alt_name=row["Alternative names"],
                cas=row["CAS"],
                smiles=row["SMILES"],
                annotation=row["Annotation"],
                comment=row["Comment"],
            )
            try:
                db.session.add(new_entry)
                db.session.commit()
            except Exception as e:
                return f"ERROR: {e}"

        return redirect(url_for("main.upload_external_from_file"))
    return render_template("upload_external_from_file.html", compounds=compounds)
