"""
User Authentication & Session Management Module:
Provides SQLite-based credential storage with SHA-256 password hashing,
user sign up, login verification, session audit logs, and user analysis history.
"""

import sqlite3
import hashlib
import os
from datetime import datetime
from typing import Optional, Dict, Any, List

DB_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "users_analytics.db")


def get_db_connection():
    os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_auth_db():
    """Initializes users and user_analyses tables."""
    conn = get_db_connection()
    cursor = conn.cursor()

    # Users table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS users (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        username TEXT UNIQUE NOT NULL,
        email TEXT UNIQUE NOT NULL,
        password_hash TEXT NOT NULL,
        full_name TEXT NOT NULL,
        institution TEXT,
        created_at TEXT NOT NULL,
        last_login TEXT
    );
    """)

    # User analyses history table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS user_analyses (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        user_id INTEGER NOT NULL,
        sample_name TEXT NOT NULL,
        disease TEXT NOT NULL,
        target_motif TEXT NOT NULL,
        seq_length INTEGER NOT NULL,
        max_repeats INTEGER NOT NULL,
        total_matches INTEGER NOT NULL,
        z_comparisons INTEGER NOT NULL,
        comparisons_saved INTEGER NOT NULL,
        prediction_label TEXT NOT NULL,
        risk_level TEXT NOT NULL,
        confidence_pct REAL NOT NULL,
        created_at TEXT NOT NULL,
        FOREIGN KEY (user_id) REFERENCES users (id)
    );
    """)

    conn.commit()
    conn.close()


def hash_password(password: str) -> str:
    """Returns SHA-256 hash of password with a static salt."""
    salt = "genomic_z_algo_salt_2026"
    return hashlib.sha256((password + salt).encode('utf-8')).hexdigest()


def register_user(username: str, email: str, password: str, full_name: str, institution: str = "") -> tuple[bool, str]:
    """Registers a new researcher/doctor user."""
    conn = get_db_connection()
    cursor = conn.cursor()

    username = username.strip().lower()
    email = email.strip().lower()

    if not username or not email or not password:
        conn.close()
        return False, "All required fields must be filled."

    pwd_hash = hash_password(password)
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    try:
        cursor.execute("""
        INSERT INTO users (username, email, password_hash, full_name, institution, created_at)
        VALUES (?, ?, ?, ?, ?, ?)
        """, (username, email, pwd_hash, full_name.strip(), institution.strip(), now_str))
        conn.commit()
        conn.close()
        return True, "Account registered successfully! You may now sign in."
    except sqlite3.IntegrityError as e:
        conn.close()
        err_msg = str(e).lower()
        if "username" in err_msg:
            return False, "Username already exists. Please choose another."
        elif "email" in err_msg:
            return False, "Email address already registered."
        return False, f"Registration error: {e}"


def authenticate_user(username: str, password: str) -> Optional[Dict[str, Any]]:
    """Authenticates username and password."""
    conn = get_db_connection()
    cursor = conn.cursor()

    username = username.strip().lower()
    pwd_hash = hash_password(password)

    cursor.execute("""
    SELECT id, username, email, full_name, institution, created_at, last_login
    FROM users
    WHERE (username = ? OR email = ?) AND password_hash = ?
    """, (username, username, pwd_hash))

    row = cursor.fetchone()
    if row:
        user_dict = dict(row)
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        cursor.execute("UPDATE users SET last_login = ? WHERE id = ?", (now_str, user_dict["id"]))
        conn.commit()
        conn.close()
        return user_dict

    conn.close()
    return None


def record_user_analysis(user_id: int, analysis_data: Dict[str, Any]):
    """Records an analysis run performed by a logged-in user."""
    conn = get_db_connection()
    cursor = conn.cursor()
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cursor.execute("""
    INSERT INTO user_analyses (
        user_id, sample_name, disease, target_motif, seq_length,
        max_repeats, total_matches, z_comparisons, comparisons_saved,
        prediction_label, risk_level, confidence_pct, created_at
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        user_id,
        analysis_data.get("sample_name", "Unknown"),
        analysis_data.get("disease", "Unknown"),
        analysis_data.get("target_motif", "CAG"),
        analysis_data.get("seq_length", 0),
        analysis_data.get("max_repeats", 0),
        analysis_data.get("total_matches", 0),
        analysis_data.get("z_comparisons", 0),
        analysis_data.get("comparisons_saved", 0),
        analysis_data.get("prediction_label", "Normal"),
        analysis_data.get("risk_level", "Low"),
        analysis_data.get("confidence_pct", 0.0),
        now_str
    ))

    conn.commit()
    conn.close()


def get_user_analyses(user_id: int) -> List[Dict[str, Any]]:
    """Retrieves all past analyses performed by a user."""
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute("""
    SELECT * FROM user_analyses
    WHERE user_id = ?
    ORDER BY id DESC
    """, (user_id,))
    rows = cursor.fetchall()
    conn.close()
    return [dict(r) for r in rows]


# Auto-initialize DB on import
init_auth_db()
