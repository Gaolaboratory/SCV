import os
import sys
import math
import time
import threading
import traceback
from multiprocessing import Process, Queue
import sqlite3
import pymol
from pymol import cmd
import zstd
from pymol2glmol import parseObjMol, parseDistObj


def create_table(db_name):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS pdbstr (
                    protein text PRIMARY KEY,
                    pdb_str blob NOT NULL,
                    ret blob NOT NULL,
                    view blob NOT NULL) ''')
    con.commit()


def is_in_db(cur, protein):
    cur.execute("SELECT * FROM pdbstr WHERE protein=?", (protein,))
    record = cur.fetchall()
    if len(record) > 0:
        return True
    return False


def insert_to_table(db_name, dict):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    try:
        cur.execute("INSERT INTO pdbstr VALUES (?, ?, ?, ?)", (dict['protein'], dict['pdb_str'], dict['ret'], dict['view']))
        con.commit()
    except sqlite3.IntegrityError:
        return


def print_progress(count, total, status=''):
    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %3.3f%% - %s" % (('='*math.ceil((count/total*50))), (count/total*100), status))
    sys.stdout.flush()


def progress_thread(q, total):
    processed = 0
    while True:
        try:
            last = q.get()
            processed += 1
            print_progress(processed, total, last)
        except Exception:
            traceback.print_exc()


def pdb2json(pdb_file):
    pdb_name = os.path.split(pdb_file)[1]
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)

    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'), os.path.join(os.getcwd(), 'js'))
        except OSError:
            pass
    try:
        cmd.set('pse_export_version', 1.74)
    except:
        pass
    names = cmd.get_session()['names']
    cmd.set('pdb_retain_ids', 1)

    ret = ''
    for obj in names:
        if (obj == None):
            continue
        if (obj[2] == 0):  # not visible
            continue
        if (obj[1] == 0 and obj[4] == 1 and obj[0] == pdb_name):
            ret += parseObjMol(obj)
            # print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    pdb_str = cmd.get_pdbstr(pdb_name)

    view_str = ''
    cmd.turn('z', 180)
    view = cmd.get_view()
    cmd.turn('z', 180)
    cx = -view[12]
    cy = -view[13]
    cz = -view[14]
    cameraZ = - view[11] - 150
    fov = float(cmd.get("field_of_view"))
    fogStart = float(cmd.get("fog_start"))
    slabNear = view[15] + view[11]
    slabFar = view[16] + view[11]
    view_str += "\nview:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % \
           (cx, cy, cz, cameraZ, slabNear, slabFar, fogStart, fov)
    for i in range(9):
        view_str += ",%.3f" % view[i]

    dict = {
            'protein': pdb_name.split('-')[1],
            'pdb_str': zstd.ZSTD_compress(pdb_str.encode('utf-8')),
            'ret': zstd.ZSTD_compress(ret.encode('utf-8')),
            'view': zstd.ZSTD_compress(view_str.encode('utf-8'))
            }
    pymol.cmd.quit()
    return dict


def batch_pdb_subroutine(db_name, path):
    sys.stdout = open(os.devnull, "w") # supress console output from pymol
    insert_to_table(db_name, pdb2json(path))


def batch_pdb_routine(q, db_name, path, pdbs):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    for index, pdb in enumerate(pdbs, start=0):
        protein = pdb.split('-')[1]
        if not is_in_db(cur, protein):
            p = Process(target=batch_pdb_subroutine, args=(db_name, os.path.join(path, pdb)))
            p.start()
            p.join()
            del p
        else:
            # print('skipping')
            pass
        try:
            q.put(pdb)
        except Queue.Full:
            traceback.print_exc()
            time.sleep(0.00001)
            q.put(pdb)
        con.commit()


def batch_pdb(db_name, path):
    pdbs = os.listdir(path)
    chunk_size = math.floor(len(pdbs)/8)
    chunks = [pdbs[i:i + chunk_size] for i in range(0, len(pdbs), chunk_size)]
    q = Queue()
    t = threading.Thread(target=progress_thread, args=(q, len(pdbs)))
    t.start()
    for chunk in chunks:
        p = Process(target=batch_pdb_routine, args=(q, db_name, path, chunk))
        p.start()


if __name__ == '__main__':
    db = sys.argv[1]
    pdb_path = sys.argv[2]

    create_table(db)
    batch_pdb(db, pdb_path)