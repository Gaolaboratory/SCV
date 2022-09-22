# main functions

import re
import pymol
import numpy as np
from pymol2glmol import *
import ahocorasick
import commons
import time


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, regex_dict=None):
    """
    generate index array from mapped peptides and PTMs

    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    from collections import defaultdict

    id_freq_array_dict = {}
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}

    regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each).replace('n', '')  # delete n in a peptide string as N term
        peptide_psm_dict[each_reg_sub].append(each)

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line_dict = {each: commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
    separtor_pos_array = commons.separator_pos(seq_line)

    aho_result = commons.automaton_matching(commons.automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1] + 1] += len(peptide_psm_dict[matched_pep])

        if ptm_index_line_dict:  # if ptm enabled
            for psm in peptide_psm_dict[matched_pep]:
                for ptm in regex_dict:
                    # print(ptm)

                    new_psm = re.sub('\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).replace('*',
                                                                                                           ptm.replace(
                                                                                                               '\\',
                                                                                                               ''))
                    # print(new_psm)
                    ptm_mod = set(re.findall(ptm, new_psm))
                    # print(ptm_mod)
                    if ptm_mod:

                        for ele in ptm_mod:
                            ### count multiple ptms in a peptide seq
                            num_of_mod = len(
                                re.findall(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm))
                            PTM_index = [m.start() for m in
                                         re.finditer(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'),
                                                     new_psm)]
                            PTM_index_clean = [ind - num * (len(ele) - 1) for ind, num in
                                               zip(PTM_index, range(num_of_mod))]
                            for indx in PTM_index_clean:
                                ptm_index_line_dict[ptm][tp[0] + indx] += 1

    for i in range(len(separtor_pos_array) - 1):
        zero_line_slice = zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
        percentage_cov = np.count_nonzero(zero_line_slice) / len(zero_line_slice) * 100
        id_freq_array_dict[id_list[i]] = zero_line_slice

        if ptm_index_line_dict:
            id_ptm_idx_dict[id_list[i]] = {ptm: np.array(
                np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[
                    0]).tolist()
                                           for ptm in ptm_index_line_dict}

    return id_freq_array_dict, id_ptm_idx_dict


def _main_(peptide_list,
           protein_dict,
           protein_id,
           pdb_file,
           ptm_color_dict=None,
           base_path=None):
    """
    output coverage 3d png or an interactive html based on GLmol

    :param peptide_list: a list of peptide with/without PTMs, e.g. ['PEPTIDE','C[142]PEPTIDE'],
    PTM should be inside a '[]'
    :param protein_dict: protein sequence dictionary, from fasta file or custom create by yourself
    :param protein_id: a uniprot ID of your interest
    :param ptm_color_dict: regex PTM and RGB color dict, e.g. {'R\[\d+\.?\d+\]':[1,1,1]}
    :param base_path: base path for html template
    :return:
    """

    time_start = time.time()
    id_freq_array_dict, id_ptm_idx_dict = freq_ptm_index_gen_batch_v2(peptide_list, protein_dict,
                                                                      regex_dict=ptm_color_dict)
    freq_arry_target = id_freq_array_dict[protein_id]
    ptm_index_target = id_ptm_idx_dict[protein_id] if id_ptm_idx_dict else None

    pdb_name = os.path.split(pdb_file)[1]
    print(pdb_name)
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)

    ## screen shot needs to be generated from pymol api
    # if pngsave:
    #     pymol.cmd.png(pngsave)
    #     print(f'image saved to {pngsave}')
    # elif not pngsave:
    #     print('image path not provided, png save failed')

    # pymol2glmol, convert pdb to pse and visualize through html
    dump_rep_color_from_array(pdb_name, freq_arry_target, ptm_index_target, ptm_color_dict, base_path)

    print(f'time used for mapping: {pdb_name} {time.time() - time_start} seconds')
    # Get out!
    pymol.cmd.quit()


if __name__ == '__main__':

    # prepare peptide list
    peptide_list = commons.modified_peptide_from_psm('D:/data/native_protein_digestion/12072021/control/0240min/psm.tsv')
    # print(peptide_list)

    # prepare protein sequence dictionary
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
    protein_dict = commons.fasta_reader(fasta_file)

    # define PTM regex color dict, avoid using red on any PTMs since default coverage is in red
    ptm_color_dict = {'n\[43\]': [0, 0, 256]}  # n terminal acetylation show in blue

    # target protein of interest
    uniprot_id = 'P61956'

    # pdb file path
    pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61956-F1-model_v1.pdb'

    # output html file path
    html_output = 'C:/Users/gao lab computer/PycharmProjects/SCV_local/P61956_test_ptm.html'

    # call main function and generate html file with interactive 3d coverage mapping
    _main_(peptide_list,
           protein_dict,
           uniprot_id,
           pdb_file=pdb_path,
           ptm_color_dict=ptm_color_dict,
           base_path=html_output)
