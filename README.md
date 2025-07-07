# COMAS Compound Submission Web Form

A Flask web application for submitting compound information for the registration to the COMAS collection. The app provides a user-friendly interface for entering user information and compound details.

## Features

- User registration with affiliation selection
- Compound data entry with validation
- Auto-assignment of plate position (MPI-DO members only)
- Submission summary and confirmation
- Data temporarily stored in SQLite databases


## Setup

1. **Clone the repository**

2. **Create a virtual environment and activate it**
   ```sh
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Install dependencies**
   ```sh
   pip install -r requirements.txt
   ```

4. **Set up environment variables**

   Ensure `.env` contains:
   ```
   FLASK-WTFS_KEY=your-secret-key
   ```

5. **Run the application**
   ```sh
   python run.py
   ```


## License

This project is for internal use.