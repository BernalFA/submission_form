# COMAS Compound Submission Web Form

A Flask web application for submitting compound information for the registration to the COMAS collection. The app provides a user-friendly interface for entering user information and compound details.

## Features

- User registration with affiliation selection
- Compound data entry with validation
- Auto-assignment of plate position (in case compounds are not served in plates)
- Submission summary and confirmation
- Data temporarily stored in SQLite databases

## Requirements

This application requires `Python 3.12`. If not already installed, please do it with

```bash
sudo apt install -y python3.12 python3.12-venv python3.12-dev python3.12-distutils
```

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

This project is for internal use only.