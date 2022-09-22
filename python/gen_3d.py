import re
import traceback

import pymol
import numpy as np
from collections import defaultdict
import sys

import zstd

from pymol2glmol import *
from glob import glob
import ahocorasick
import commons
import sys
import json
import sqlite3
import os


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
        assert seq_line[start_idx:start_idx + len(original_value)] == original_value
    return result


def fasta_reader(fasta_file_path):
    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def peptide_counting(peptide_tsv_file):
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string, regex_pat='\w{1}\[\d+\.?\d+\]'):
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

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list, PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:

            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) - 1
            # print (start_pos,end_pos,new_pep)
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:
                    PTM_index = pep.find(ele)
                    # PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos + PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)

    return freq_array, PTM_loc_list, PTM_sites_counting


def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]', line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list


def color_getter(freq_array, ptm_idx_dict, regex_color_dict, pdb_str):
    """
    based on numpy array get blocks of zeros and non_zeros, generate color string for GLMOL
    :param freq_array: numpy array
    :param pdb_str: cmd.get_pdbstr()
    :return:
    """
    import re
    import numpy as np
    from collections import defaultdict

    # from pdb_str get element pos to amino acid pos
    pdb_str = pdb_str.split('\nTER')[0].split('\n')
    amino_ele_pos_dict = defaultdict(list)
    for line in pdb_str:
        amino_ele_pos_dict[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append(
            int(re.search('\d+', line).group()))
    amino_ele_pos_dict = {each: sorted(amino_ele_pos_dict[each]) for each in amino_ele_pos_dict}

    defalt = 'color:0.500,0.500,0.500:'  # grey
    covered = 'color:1.000,0.000,0.000:'
    # get index of zeros and nonzeros
    non_zero_index = np.nonzero(freq_array)[0]
    zero_index = np.nonzero(freq_array == 0)[0]

    # get index blocks of zeros and nonzeros
    cov_pos_block = np.split(non_zero_index, np.where(np.diff(non_zero_index) != 1)[0] + 1)
    non_cov_pos_block = np.split(zero_index, np.where(np.diff(zero_index) != 1)[0] + 1)
    print(cov_pos_block, non_cov_pos_block)

    # string concatenate
    defalt += ','.join([str(amino_ele_pos_dict[each[0] + 1][0]) + '-' + str(amino_ele_pos_dict[each[-1] + 1][-1])
                        for each in non_cov_pos_block])
    covered += ','.join([str(amino_ele_pos_dict[each[0] + 1][0]) + '-' + str(amino_ele_pos_dict[each[-1] + 1][-1])
                         for each in cov_pos_block])

    # ptm color string concatenate
    ptm_color = ''
    if ptm_idx_dict:
        for ptm in ptm_idx_dict:
            ptm_color += 'color:' + ','.join(
                ['%.3f' % (int(each) / 256) for each in regex_color_dict[ptm].strip('][').split(', ')]) + ':'
            ptm_color += ','.join([str(amino_ele_pos_dict[idx + 1][0]) + '-'
                                   + str(amino_ele_pos_dict[idx + 1][-1])
                                   for idx in ptm_idx_dict[ptm]])
            ptm_color += '\n'
    return defalt + '\n' + covered + '\n' + ptm_color.rstrip('\n')


def dump_rep_server(name, freq_array, ptm_idx_dict, regex_color_dict, bg_color):
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'),
                            os.path.join(os.getcwd(), 'js'))
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
        if (obj[1] == 0 and obj[4] == 1 and obj[0] == name):
            ret += parseObjMol(obj)
            print(ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    start = time.time()
    pdb_str = cmd.get_pdbstr(name)
    print(f'time used for pbdstr: {time.time() - start}')

    start = time.time()
    ret += '\n' + color_getter(freq_array, ptm_idx_dict, regex_color_dict, pdb_str)
    print(f'time used for color_getter: {time.time() - start}')

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
    ret += "\nview:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % \
           (cx, cy, cz, cameraZ, slabNear, slabFar, fogStart, fov)
    for i in range(9):
        ret += ",%.3f" % view[i]

    ret += "\nbgcolor:" + bg_color
    # bgcolor = cmd.get_setting_tuple('bg_rgb')[1]
    #
    # if len(bgcolor) == 1:
    #     bgcolor = cmd.get_color_tuple(bgcolor[0])
    #
    # ret += "\nbgcolor:%02x%02x%02x" % (int(255 * float(bgcolor[0])), \
    #                                    int(255 * float(bgcolor[1])), int(255 * float(bgcolor[2])))
    # if 'PYMOL_GIT_MOD' in os.environ:
    #     template = open(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'imported.html')).read().\
    #         replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
    #         replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    # else:
    #     template = open('imported.html').read().\
    #         replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
    #         replace('###INCLUDE_REPRESENTATION_HERE###',ret)
    #
    # if base_path:
    #     f = open(base_path+name + '.html', 'w')
    #     print (f'html file to {base_path+name}.html')
    # else:
    #     f = open(name.split('-')[1] + '.html', 'w')
    #     print ('html file to %s' % name.split('-')[1]+'.html')
    # f.write(template)
    # f.close()

    #     dict = {'pbdstr': cmd.get_pdbstr(name), 'ret': ret}
    dict = {'pbdstr': pdb_str, 'ret': ret}

    # with open(base_path + '.json', 'w') as f:
    #     json.dump(dict, f)
    # pymol.cmd.quit()

    return dict


def show_cov_3d_v2(protein_id,
                   job_number,
                   pdb_file,
                   frequency_array,
                   id_ptm_idx_dict,
                   regex_color_dict=None,
                   png_sava_path=None,
                   base_path=None,
                   background_color='black'):
    """

    :param peptide_list:
    :param protein_seq:
    :param pdb_file:
    :param id_freq_array_dict: returned by freq_ptm_index_gen_batch_v2
    :param id_ptm_idx_dict: returned by freq_ptm_index_gen_batch_v2
    :param png_sava_path:
    :param base_path: html output base path
    :return:
    """
    time_start = time.time()
    if id_ptm_idx_dict != {}:
        ptm_nonzero_idx_dict = id_ptm_idx_dict
    else:
        ptm_nonzero_idx_dict = None

    pdb_name = os.path.split(pdb_file)[1]
    print(pdb_name)
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable()
    print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('black')
    print(f'time used for preproccessing: {time.time() - time_start}')

    # set customized color
    #     if regex_color_dict:
    #         for i in regex_color_dict:
    #             print(regex_color_dict[i])
    #             pymol.cmd.set_color(i, regex_color_dict[i])
    #
    #     print (ptm_nonzero_idx_dict)
    #     max_freq = np.max(frequency_array)
    #     for i in range(len(frequency_array)): # iterate over each residue position
    #         ptm = False
    #         if ptm_nonzero_idx_dict:
    #             for ptm_regex in ptm_nonzero_idx_dict:
    #                 # if index has ptm
    #                 if i in ptm_nonzero_idx_dict[ptm_regex]:
    #                     ptm = True
    #                     pymol.cmd.color(ptm_regex, 'resi %i' % (i + 1))
    #
    #         if ptm == False:  # if no ptm found
    #             if frequency_array[i] == 0:
    #                 pymol.cmd.color('grey', 'resi %i' % (i + 1))
    #            # elif 1 <= frequency_array[i] < 0.2 * max_freq:
    #            #     pymol.cmd.color('paleyellow', 'resi %i' % (i + 1))
    #            # elif 0.2 * max_freq <= frequency_array[i] < 0.4 * max_freq:
    #            #     pymol.cmd.color('tv_yellow', 'resi %i' % (i + 1))
    #            # elif 0.4 * max_freq <= frequency_array[i] < 0.6 * max_freq:
    #            #     pymol.cmd.color('yelloworange', 'resi %i' % (i + 1))
    #            # elif 0.6 * max_freq <= frequency_array[i] < 0.8 * max_freq:
    #            #     pymol.cmd.color('tv_orange', 'resi %i' % (i + 1))
    #            # else:
    #            #    pymol.cmd.color('sand', 'resi %i' % (i + 1))
    #             else:
    #                 pymol.cmd.color('red','resi %i' % (i + 1))
    #         else: # ptm color assigned, move on to the next residue
    #             continue

    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print(f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    #     dump_rep(pdb_name,base_path)
    #     dump_rep_2(pdb_name, os.path.join(os.getcwd(), 'results/' + job_number + '_' + protein_id), background_color)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # Get out!
    # pymol.cmd.quit()
    return dump_rep_server(pdb_name, frequency_array, id_ptm_idx_dict, regex_color_dict, background_color)


def pdb_file_reader(pdb_file):
    """
    reads a pdb file into protein sequence
    :param pdb_file:
    :return:
    """
    aa_dict = {'ALA': 'A',
               'ARG': 'R',
               'ASN': 'N',
               'ASP': 'D',
               'ASX': 'B',
               'CYS': 'C',
               'GLU': 'E',
               'GLN': 'Q',
               'GLX': 'Z',
               'GLY': 'G',
               'HIS': 'H',
               'ILE': 'I',
               'LEU': 'L',
               'LYS': 'K',
               'MET': 'M',
               'PHE': 'F',
               'PRO': 'P',
               'SER': 'S',
               'THR': 'T',
               'TRP': 'W',
               'TYR': 'Y',
               'VAL': 'V'}

    aa_reg_str = '|'.join([key for key in aa_dict])

    import re

    with open(pdb_file, 'r') as f_o:
        f_split = f_o.read().split('\nATOM')[1:]
        aa_list = [re.search(aa_reg_str, each).group(0) for each in f_split]
        aa_list = [aa_dict[each] for each in aa_list]
    protein_seq = ''
    for i in range(len(aa_list) - 1):
        if aa_list[i + 1] == aa_list[i]:
            continue
        else:
            protein_seq += aa_list[i]
    # add last aa
    protein_seq += aa_list[-1]
    return protein_seq


def is_in_db(cur, protein):
    cur.execute("SELECT * FROM pdbstr WHERE protein=?", (protein,))
    record = cur.fetchall()
    if len(record) > 0:
        return True
    return False


def get_from_db(cur, protein):
    cur.execute("SELECT * FROM pdbstr WHERE protein=?", (protein,))
    record = cur.fetchall()
    return {
        'protein': record[0][0],
        'pdb_str': zstd.ZSTD_uncompress(record[0][1]).decode('utf-8'),
        'ret': zstd.ZSTD_uncompress(record[0][2]).decode('utf-8'),
        'view': zstd.ZSTD_uncompress(record[0][3]).decode('utf-8')
    }


if __name__ == '__main__':
    import time
    import traceback
    import pandas as pd
    import zmq

    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind('ipc:///tmp/gen3d.ipc')
    print('bound on ipc:///tmp/gen3d.ipc')

    while True:
        message = socket.recv()
        start = time.time()
        try:
            recv_dict = json.loads(message)
            if os.path.splitext(recv_dict['pdb_dest'])[1] == '.db':
                con = sqlite3.connect(recv_dict['pdb_dest'])
                cur = con.cursor()
                if is_in_db(cur, recv_dict['protein']):
                    db_res = get_from_db(cur, recv_dict['protein'])
                    print(recv_dict['freqArr'])
                    result_dict = {
                        'job_number': recv_dict['job_number'],
                        'pdbstr': db_res['pdb_str'],
                        'ret': db_res['ret'] + "\n" +
                               color_getter(
                                   np.array(recv_dict['freqArr'].split(',')).astype(np.int32),
                                   json.loads(recv_dict['id_ptm_idx_dict']),
                                   json.loads(recv_dict['regex_dict']), db_res['pdb_str'])
                               + db_res['view']
                               + "\nbgcolor:" + recv_dict['background_color'] + "\n"
                    }
                    print(f'time used for processing: {time.time() - start}')
                    socket.send_string(json.dumps(result_dict))
                else:
                    socket.send_string("error")
                con.close()
            else:
                for pdb in os.listdir(recv_dict['pdb_dest']):
                    if pdb == recv_dict['protein']:
                        dump_dict = show_cov_3d_v2(recv_dict['protein'], recv_dict['job_number'],
                                                   os.path.join(recv_dict['pdb_dest'], pdb),
                                                   np.array(recv_dict['freqArr'].split(',')).astype(np.int32),
                                                   id_ptm_idx_dict=json.loads(recv_dict['id_ptm_idx_dict']),
                                                   regex_color_dict=json.loads(recv_dict['regex_dict']),
                                                   background_color=recv_dict['background_color'])
                        result_dict = {
                            'job_number': recv_dict['job_number'],
                            'pdbstr': dump_dict['pbdstr'],
                            'ret': dump_dict['ret']
                        }
                        socket.send_string(json.dumps(result_dict))
        except Exception:
            traceback.print_exc()
            socket.send_string("error")
        # socket.send_string(str(is_in_db(cur, message.decode('utf-8'))))

#     job_number = sys.argv[1]
#     print(job_number)
#
#     percent = float(sys.argv[2])
#     print(percent)
#     protein = sys.argv[3]
#     print(protein)
#     freqArr = np.array(sys.argv[4].split(',')).astype(np.int32)
#     print(freqArr)
#     print(sys.argv[5])
#     id_ptm_idx_dict = json.loads(sys.argv[5])
#     print(sys.argv[6])
#     regex_dict = json.loads(sys.argv[6])
#     background_color = sys.argv[7]
#     pdb_dest = sys.argv[8]
# if os.path.exists(pdb_dest + '/AF-'+ protein + '-F1-model_v1.pdb'):
# #if os.path.exists(pdb_dest + '/AF-'+ protein.split('_')[1] + '-F1-model_v1.pdb'):
#     print('exists!')
#     show_cov_3d_v2(protein, job_number, pdb_dest + '/AF-'+ protein + '-F1-model_v1.pdb',
#                             freqArr, id_ptm_idx_dict=id_ptm_idx_dict,regex_color_dict=regex_dict, background_color=background_color)
# else:
#     print('doesnt exists')
#     for pdb in os.listdir(pdb_dest):
#         if pdb == protein:
#            show_cov_3d_v2(protein, job_number, os.path.join(pdb_dest, pdb),
#                             freqArr, id_ptm_idx_dict=id_ptm_idx_dict,regex_color_dict=regex_dict, background_color=background_color)
#            break


# if is_in_db(cur, protein):
#     start = time.time()
#     dict = get_from_db(cur, protein)
#     print(f'time used for get_from_db: {time.time() - start}')
#     color_getter(freqArr, id_ptm_idx_dict, regex_dict, dict['pdb_str'])
#     print(f'time used for color_getter: {time.time() - start}')
#
# else:
#     print('not yet implemented')
