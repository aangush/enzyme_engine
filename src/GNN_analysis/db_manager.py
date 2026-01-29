import sqlite3
import os

class ScoutDB:
    def __init__(self, db_path="data/scout.db"):
        self.db_path = db_path
        self._initialize_db()

    def _initialize_db(self):
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE IF NOT EXISTS anchors (protein_id TEXT PRIMARY KEY, product TEXT, organism TEXT, sequence TEXT)''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS instances (id INTEGER PRIMARY KEY AUTOINCREMENT, anchor_id TEXT, nuc_accession TEXT, start_pos INTEGER, end_pos INTEGER, strand INTEGER, UNIQUE(anchor_id, nuc_accession, start_pos))''')
        cursor.execute('''CREATE TABLE IF NOT EXISTS neighbors (id INTEGER PRIMARY KEY AUTOINCREMENT, instance_id INTEGER, protein_id TEXT, product TEXT, distance_bp INTEGER, direction TEXT, sequence TEXT, FOREIGN KEY(instance_id) REFERENCES instances(id))''')
        conn.commit()
        conn.close()

    def instance_exists(self, anchor_id, nuc_acc):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT id FROM instances WHERE anchor_id = ? AND nuc_accession = ?", (anchor_id, nuc_acc))
        result = cursor.fetchone()
        conn.close()
        return result[0] if result else None

    def add_instance(self, anchor_id, nuc_acc, start, end, strand):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("INSERT OR IGNORE INTO instances (anchor_id, nuc_accession, start_pos, end_pos, strand) VALUES (?, ?, ?, ?, ?)",
                       (anchor_id, nuc_acc, start, end, strand))
        last_id = cursor.lastrowid
        if last_id == 0: # If IGNORE triggered, find existing ID
            cursor.execute("SELECT id FROM instances WHERE anchor_id = ? AND nuc_accession = ?", (anchor_id, nuc_acc))
            last_id = cursor.fetchone()[0]
        conn.commit()
        conn.close()
        return last_id

    def add_neighbor(self, instance_id, pid, product, dist, direction, seq=""):
        conn = sqlite3.connect(self.db_path)
        conn.execute("INSERT INTO neighbors (instance_id, protein_id, product, distance_bp, direction, sequence) VALUES (?, ?, ?, ?, ?, ?)",
                     (instance_id, pid, product, dist, direction, seq))
        conn.commit()
        conn.close()