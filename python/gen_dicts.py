import re
import traceback

import pymol
import numpy as np
from collections import defaultdict
import sys
from pymol2glmol import *
from glob import glob
import ahocorasick
import commons
import sys
import json
import time
import sqlite3
import zmq


pdb_dict = {
            'mouse': '/home/cgrams/Protein-Vis/db/UP000000589_10090_MOUSE.db',
            'rat': '/home/cgrams/Protein-Vis/db/UP000002494_10116_RAT.db',
            'human': '/home/cgrams/Protein-Vis/db/UP000005640_9606_HUMAN.db'
            }

fasta_dict = {
            'mouse': '/home/cgrams/Protein-Vis/fastas/uniprot-proteome_UP000000589_sp_only_mouse.fasta',
            'rat': '/home/cgrams/Protein-Vis/fastas/uniprot-proteome_UP000002494_sp_only_rat.fasta',
            'human': '/home/cgrams/Protein-Vis/fastas/uniprot-proteome_UP000005640_sp_only_human.fasta'
            }


def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A


def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx+len(original_value)] == original_value
    return result


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def fasta_reader2(fasta_file_path):
   # use gene_uniprot as key
    gene_protein_seq_dict = {}
    with open(fasta_file_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        gene_protein_seq_dict[gene+'_'+uniprot_id] = seq
    return gene_protein_seq_dict


def fasta_reader3(fasta_path:str):

    protein_dict = {}
    with open(fasta_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        des = ' '.join(first_line.split(' ')[1:]).split(' OS=')[0]
        protein_dict[uniprot_id] = (seq,gene,des)
    return protein_dict

def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string,regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map single protein seq
    :param peptide_list:
    :param protein_seq_string:
    :param regex_pat:
    :return:
    """

    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:

            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            # print (start_pos,end_pos,new_pep)
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:

                    PTM_index = pep.find(ele)
                   #PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos+PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)

    return freq_array, PTM_loc_list, PTM_sites_counting


def freq_ptm_index_gen_batch(psm_list, protein_dict, regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map large psm list on whole proteome
    :param psm_list: psm list with ptms
    :param protein_dict: protein sequence dictionary
    :param regex_pat: regex pattern
    :return:
    """
    from collections import Counter
    from collections import defaultdict
    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}

    peptide_psm_dict = defaultdict(list) # append all psm into a dictionary
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line = commons.zero_line_for_seq(seq_line)
    separtor_pos_array = commons.separator_pos(seq_line)

    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]),seq_line)

    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
        # print (tp[0],tp[1],matched_pep, 'aho')
        for psm in peptide_psm_dict[matched_pep]:
            psm_mod = re.findall(regex_pat,psm)
            if psm_mod: # if psm has mod
                for ele in psm_mod:
                    ptm_idx = psm.find(ele)
                    ptm_index_line[tp[0]+ptm_idx]+=1

    for i in range(len(separtor_pos_array)-1):

        id_freq_array_dict[id_list[i]] = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
        id_ptm_idx_dict[id_list[i]] = np.nonzero(ptm_index_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1

    return id_freq_array_dict,id_ptm_idx_dict


# def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, regex_dict=None):
#     """
#
#     :param psm_list:
#     :param protein_dict:
#     :param regex_dict: {regex:HEX color}
#     :return:
#     """
#
#     from collections import Counter
#     from collections import defaultdict
#     id_freq_array_dict = {}
#     ptm_site_counting = defaultdict(int)
#     id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
#
#     regex_pat = '\w{1}\[\d+\.?\d+\]' # universal ptm pattern
#     peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
#     for each in psm_list:
#
#         each_reg_sub = re.sub(regex_pat, my_replace, each)
#         peptide_psm_dict[each_reg_sub].append(each)
#
#     # aho mapping
#     id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
#     seq_line = commons.creat_total_seq_line(seq_list, sep="|")
#     zero_line = commons.zero_line_for_seq(seq_line)
#     ptm_index_line_dict = {each:commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
#     separtor_pos_array = commons.separator_pos(seq_line)
#
#     aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
#     for tp in aho_result:
#         matched_pep = tp[2]  # without ptm site
#         zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
#         # ptm assign might need to be optimized
#         if ptm_index_line_dict:  # if ptm enabled
#             for psm in peptide_psm_dict[matched_pep]:
#                 for ptm in regex_dict:
#                     # only keep one ptm in psm if there are multiple for correct index finding
#                     new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)
#
#                     ptm_mod = re.findall(ptm, new_psm)
#                     print(ptm_mod)
#                     if ptm_mod:
#                         for ele in ptm_mod:
#                             print (new_psm)
#                             ptm_idx = new_psm.find(ele)
#                             print(matched_pep, tp[0], ele, ptm_idx)
#                             ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1
#
#     for i in range(len(separtor_pos_array)-1):
#         zl = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
#         # id_freq_array_dict[id_list[i]] = []
#         for elem in zl:
#             if elem != 0:
#                 id_freq_array_dict[id_list[i]] = zl
#                 break
#
#         # id_freq_array_dict.pop(id_list[i])
#         # id_freq_array_dict[id_list[i]] = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
#         if ptm_index_line_dict:
#             id_ptm_idx_dict[id_list[i]]= {ptm:np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1
#                                           for ptm in ptm_index_line_dict}
#     print (id_ptm_idx_dict)
#     return id_freq_array_dict, id_ptm_idx_dict

def pdb_file_reader(pdb_file_list:list):
    """
    reads a pdb file into protein sequence
    :param pdb_file:
    :return:
    """
    aa_dict = {'ALA':'A',
               'ARG':'R',
               'ASN':'N',
               'ASP':'D',
               'ASX':'B',
               'CYS':'C',
               'GLU':'E',
               'GLN':'Q',
               'GLX':'Z',
               'GLY':'G',
               'HIS':'H',
               'ILE':'I',
               'LEU':'L',
               'LYS':'K',
               'MET':'M',
               'PHE':'F',
               'PRO':'P',
               'SER':'S',
               'THR':'T',
               'TRP':'W',
               'TYR':'Y',
               'VAL':'V'}

    aa_reg_str = '|'.join([key for key in aa_dict])

    import re
    import os


    pdb_protein_seq_dict = {}

    for pdb_file in pdb_file_list:
        with open(pdb_file,'r') as f_o:
            f_split = f_o.read().split('\nATOM')[1:]
            pos_aa_list = [(int(re.search('\d+(?=\s+[+-]?\d+\.)',each).group()),
                           re.search(aa_reg_str,each).group(0)) for each in f_split]
            protein_seq = ''
            for i in range(len(pos_aa_list)-1):
                if pos_aa_list[i+1][0] == pos_aa_list[i][0]:
                    continue
                else:
                    protein_seq += aa_dict[pos_aa_list[i][1]]
    
        # add last aa
        protein_seq+=aa_dict[pos_aa_list[-1][1]]
        pdb_protein_seq_dict[os.path.split(pdb_file)[-1]] = (protein_seq, )
        print (protein_seq)
    return pdb_protein_seq_dict


def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, pdbs, regex_dict=None):
    """
    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """
    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop


    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    h = []
    regex_pat = '\w{1}\[\d+\.?\d+\]' # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)
    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line_dict = {each:commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
    separtor_pos_array = commons.separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
        # ptm assign might need to be optimized
        if ptm_index_line_dict:  # if ptm enabled
            print(regex_dict)
            for psm in peptide_psm_dict[matched_pep]:
                for ptm in regex_dict:
#                     print(ptm)
                    # only keep one ptm in psm if there are multiple for correct index finding
                    new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)
                    ptm_mod = re.findall(ptm, new_psm)
#                     print(ptm_mod)
                    if ptm_mod:
#                        for ele in ptm_mod:
#                            print (new_psm)
#                            ptm_idx = new_psm.find(ele)
#                            print(matched_pep, tp[0], ele, ptm_idx)
#                            ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1
                        for ele in ptm_mod:
                            num_of_mod = len(re.findall(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm))
                            PTM_index = [m.start() for m in re.finditer(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm)]
                            PTM_index_clean = [ind - num * (len(ele) - 1) for ind, num in zip(PTM_index, range(num_of_mod))]
                            for indx in PTM_index_clean:
                                ptm_index_line_dict[ptm][tp[0] + indx] += 1
    time_start = time.time()
    for i in range(len(separtor_pos_array)-1):
        zero_line_slice = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
        if len(zero_line_slice) > 0:
            percentage_cov = np.count_nonzero(zero_line_slice)/len(zero_line_slice)*100
        else:
            percentage_cov = 0.0
        if percentage_cov != 0.0:
            heappush(h,(percentage_cov,(id_list[i], protein_dict[id_list[i]][1:]),zero_line_slice.tolist(), any(id_list[i] in pdb for pdb in pdbs)))
            id_freq_array_dict[id_list[i]] = zero_line_slice.tolist()

        if ptm_index_line_dict:
            id_ptm_idx_dict[id_list[i]]= {ptm:np.array(np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]).tolist()
                                                                                   for ptm in ptm_index_line_dict}
    print (time.time()-time_start)

    return id_freq_array_dict, id_ptm_idx_dict, [heappop(h) for i in range(len(h))][::-1]


def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]',line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list


def fetch_all_proteins_from_db(cur):
    cur.execute("SELECT protein FROM pdbstr")
    records = cur.fetchall()
    return records


def create_table(db_name):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS results (
                    job_number text PRIMARY KEY,
                    pq text NOT NULL,
                    id_ptm_idx_dict text NOT NULL,
                    regex_dict text NOT NULL,
                    background_color text NOT NULL,
                    pdb_dest text NOT NULL) ''')
    con.commit()


def insert_to_table(db_name, job_number, pq, id_ptm_idx_dict, regex_dict, background_color, pdb_dest):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute("INSERT INTO results VALUES (?, ?, ?, ?, ?, ?)",
                (job_number,
                 json.dumps(pq),
                 json.dumps(id_ptm_idx_dict),
                 json.dumps(regex_dict),
                 background_color,
                 pdb_dest)
                )
    con.commit()


if __name__=='__main__': 
    import pandas as pd
#     print(sys.argv)
    result_db = '../db/SCV.db'
    create_table(result_db)
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind('ipc:///tmp/genDicts.ipc')
    print('bound on ipc:///tmp/genDicts.ipc')

    while True:
        message = socket.recv()
        try:
            recv_dict = json.loads(message)

            job_number = recv_dict['job_number']
            print(job_number)
            PSMs = recv_dict['psms'].split(',')
            PTMs = {}

            if recv_dict['ptms'] != '':
                PTMs = json.loads(recv_dict['ptms'])

            regex_dict = {}
            for key in PTMs:
                regex_dict[key.replace('[', '\[').replace(']','\]')] = str(PTMs[key])

            background_color = recv_dict['background_color']
            pdb_dest = recv_dict['pdb_dest']
            species = recv_dict['species']
            start = time.time()
            pdbs = set()
            if pdb_dest == '':
                pdb_dest = pdb_dict[species]
                con = sqlite3.connect(pdb_dest)
                con.row_factory = lambda cursor, row: row[0]
                cur = con.cursor()
                protein_dict = fasta_reader3(fasta_dict[species])
                db_proteins = set(fetch_all_proteins_from_db(cur))
                for protein in protein_dict:
                    if protein in db_proteins:
                        pdbs.add(protein)
            else:
                pdb_dest = os.path.join('/home/cgrams/Protein-Vis/', pdb_dest+'/')
                files = os.listdir(pdb_dest)
                for file in files:
                    pdbs.add((os.path.join(pdb_dest, file)))
                protein_dict = pdb_file_reader(pdbs)
            print(PSMs)
            print(protein_dict)
            print(pdbs)
            print(regex_dict)
            id_freq_array_dict, id_ptm_idx_dict, pq = freq_ptm_index_gen_batch_v2(PSMs,protein_dict, pdbs, regex_dict=regex_dict)

            print('coverage: ' + str(len(id_freq_array_dict)))

            insert_to_table(result_db, job_number, pq, id_ptm_idx_dict, regex_dict, background_color, pdb_dest)

            socket.send(b'ok')
        except Exception:
            traceback.print_exc()
            socket.send(b'error')
