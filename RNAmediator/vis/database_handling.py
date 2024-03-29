import base64
import csv
import gzip
import io
import os
import sqlite3
import time
from itertools import groupby
from operator import itemgetter
from typing import List

import numpy as np
import pandas as pd

COLUMN_NAMES = {
    'Chr': "VARCHAR(10)",
    'Start': "INT",
    'End': "INT",
    'Accessibility_difference': "REAL",
    'Strand': "VARCHAR(2)",
    'Distance_to_constraint': "INT",
    'Accessibility_no_constraint': "REAL",
    'Accessibility_constraint': "REAL",
    'Energy_Difference': "REAL",
    'Kd_change': "REAL",
    'Zscore': "REAL",
    'Genomic_Start': 'INT',
    'Genomic_END': 'INT',
    'Gene_of_interest': 'VARCHAR(20)',
}


class SearchSettings:
    def __init__(self, chromosome="", goi="", span_start=0,
                 span_end="Infinity"):
        self.chr = chromosome if chromosome is not None else ""
        self.goi = goi if goi is not None else ""
        try:
            self.span_start = int(span_start)
        except (ValueError, TypeError):
            self.span_start = 0
        try:
            self.span_end = int(span_end)
        except (ValueError, TypeError):
            self.span_end = "Infinity"

    def __str__(self):
        return str(self.__dict__)


def get_interesting(
        db_path: str,
        page: int = 0,
        ordering: str = "Max_Value",
        sorting_clicks: int = 0,
        substrings: SearchSettings = None,
        number_of_interesting: int = 10,
        intersect_filter: List[str] = None
):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    if sorting_clicks % 2:
        deasc = "ASC"
    else:
        deasc = "DESC"
    if substrings is None:
        substrings = SearchSettings()
    cur = conn.cursor()
    return_list = []
    if intersect_filter is not None and "ALL" not in intersect_filter and len(intersect_filter) > 0:
        for entry in intersect_filter:
            cur.execute(
                f"SELECT * FROM importance "
                f"WHERE Chr like '{substrings.chr}%' "
                f"AND Gene_of_interest like '{substrings.goi}%' "
                f"AND Genomic_Start >= ? "
                f"AND Genomic_End <= ? "
                f"AND Gene_of_interest IN (SELECT Gene_of_Interest FROM {entry})"
                f"ORDER BY {ordering} {deasc} "
                f"LIMIT {number_of_interesting} OFFSET {number_of_interesting * page}",
                (substrings.span_start, substrings.span_end),
            )
            return_list += cur.fetchall()
    else:
        cur.execute(
            f"SELECT * FROM importance "
            f"WHERE Chr like '{substrings.chr}%' "
            f"AND Gene_of_interest like '{substrings.goi}%' "
            f"AND Genomic_Start >= ? "
            f"AND Genomic_End <= ? "
            f"ORDER BY {ordering} {deasc} "
            f"LIMIT {number_of_interesting} OFFSET {number_of_interesting * page}",
            (substrings.span_start, substrings.span_end),
        )
        return_list = cur.fetchall()
    conn.close()
    return return_list


def csv_to_sqlite(file: str, db_path: str):
    assert not os.path.exists(db_path), "database already exists."
    con = sqlite3.connect(db_path)
    con.execute('PRAGMA synchronous = OFF')
    cur = con.cursor()
    names = [f"{key} {COLUMN_NAMES[key]}" for key in COLUMN_NAMES]
    columns = ", ".join(names)
    call = f'CREATE TABLE IF NOT EXISTS test ( {columns} ) '
    cur.execute(call)
    if ".gz" in file:
        file_handle = gzip.open(file, "rt")
    else:
        file_handle = open(file)
    cur.executemany('INSERT INTO test VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);',
                    insert_generator(file_handle))
    file_handle.close()
    cur.execute(
        "CREATE INDEX test_idx ON test (Chr, Gene_of_interest, Genomic_Start, Genomic_End)")
    file_handle.close()
    con.commit()
    con.close()
    insert_interesting_table(db_path)


def insert_intersect(contents, filename, date, db_path):
    conn = sqlite3.connect(db_path)
    table = filename.split(".")[0]
    cur = conn.cursor()
    cur.execute(
        ''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name=? ''',
        [table])
    if not cur.fetchone()[0] == 1:
        cur.execute(f'CREATE TABLE IF NOT EXISTS {table} '
                    f'( Chr VARCHAR(5), '
                    f'Genomic_Start INT, '
                    f'Genomic_End INT, '
                    f'Gene_of_Interest VARCHAR(20), '
                    f'Distance INT, '
                    f'Strand VARCHAR(2)) ')
    try:
        cur.executemany(
            f'INSERT INTO {table} VALUES (?,?,?,?,?,?);',
            intersect_generator(contents))
    except (AssertionError, ValueError):
        cur.execute(
            f"DROP TABLE IF EXISTS {table};"
        )
        raise(sqlite3.DatabaseError)
    finally:
        conn.commit()
        conn.close()


def intersect_generator(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    assert len(decoded) > 0
    datastream = io.StringIO(decoded.decode('utf-8'))
    csv_reader = csv.reader(
        datastream,
        delimiter="\t",
    )
    for row in csv_reader:
        chrom, start, stop, gene_name, _, strand, distance = row[0:7]
        gene_name = gene_name.split("|")[0]
        assert strand in ["+", "-"]
        assert int(distance)

        yield [chrom, start, stop, gene_name, distance, strand]


def get_intersects(database, genomic_start, genomic_end, strand):
    con = sqlite3.connect(database)
    cursor = con.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    names = [name[0] for name in cursor.fetchall()]
    names.remove("test")
    names.remove("importance")
    tuples = []
    for name in names:
        cursor.execute(
            f"SELECT Distance FROM {name} "
            "WHERE Genomic_Start >= ? "
            "AND Genomic_End <=? "
            "AND Strand = ?",
            (genomic_start, genomic_end, strand)
        )
        distances = list(set([x[0] for x in cursor.fetchall()]))
        distances.sort()
        output = []
        for k, g in groupby(enumerate(distances), lambda ix: ix[0] - ix[1]):
            output.append([*map(itemgetter(1), g)])
        tuples.append((name, output))
    return tuples


def insert_generator(file_handle):
    csv_reader = csv.reader(
        file_handle,
        delimiter="\t",
    )
    for x, row in enumerate(csv_reader):
        goi, pos, genomic_pos = row[3].split("|")
        gstart, gend = genomic_pos.split("-")
        row = row[:-1]
        row = row + [int(gstart), int(gend), goi]
        row.pop(3)
        yield row


def insert_interesting_table(db_path: str):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute(f"CREATE TABLE IF NOT EXISTS importance "
                f"( Chr VARCHAR(5), "
                f"Gene_of_interest VARCHAR(20), "
                f"Genomic_Start INT,"
                f"Genomic_End INT, "
                f"Mean_Value REAL, "
                f"Max_Value REAL, "
                f"Max_Zscore REAL, "
                f" FOREIGN KEY (Gene_of_interest) REFERENCES test(Gene_of_interest)) ")

    cur.execute(
        "SELECT DISTINCT Chr, Gene_of_interest, Genomic_Start, Genomic_End FROM test")
    constraints = cur.fetchall()

    for entry in constraints:
        cur.execute(
            "SELECT Distance_to_constraint, Accessibility_difference, Chr, Zscore FROM test "
            "WHERE Chr=? AND Gene_of_interest=? AND Genomic_Start=? AND Genomic_End=?",
            entry)
        values = cur.fetchall()
        constraint_max = np.max([abs(diff[1]) for diff in values])
        constraint_mean = np.mean([abs(diff[1]) for diff in values])
        zscore_max = np.max([abs(score[3]) for score in values])
        cur.execute("INSERT INTO importance VALUES (?, ?, ?, ?, ?, ?, ?)",
                    [entry[0], entry[1], entry[2], entry[3],
                     constraint_mean, constraint_max, zscore_max])
    cur.execute("CREATE INDEX interesting_idx "
                "ON importance (Chr, Gene_of_interest, Genomic_Start, Genomic_End, Max_Value, Mean_Value, Max_Zscore)")
    con.commit()
    con.close()


def read_data(file: str):
    return pd.read_csv(
        file, delimiter='\t', names=list(COLUMN_NAMES)
    )  # ,'ChrBS','StartBS','EndBS','NameBS','ScoreBS','StrandBS'])
