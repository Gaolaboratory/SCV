![logo](/images/scv_logo.png)

## Sequence coverage visualizer (SCV)
Generate an interactive 3D protein coverage map (PTM enabled) given peptides list and protein sequence

Developed and maintained by Xinhao shao and Chris Grams at Gao lab, UIC. For any bugs or issues, please contact xshao8@uic.edu

### Usage

#### SCV web application:
[scv.lab.gy](http://scv.lab.gy/)

tutorial: [https://youtu.be/NVoVMHVczKY](https://youtu.be/NVoVMHVczKY)

Example using uploaded .pdb file from [RCSB](https://www.rcsb.org/) other than default Alphafold pdb structures:

1. Download a monomeric structure [3hff](https://www.rcsb.org/structure/3HFF) from RCSB database 
2. Upload downloaded 3hff.pdb in [SCV](http://scv.lab.gy/)
3. Copy and paste PSMs of your interest (e.g. ATKAVAVLKGDGPVQG)
4. Click 'Next' and visualize the sequence coverage on your own PDB structure with [SCV](http://scv.lab.gy/)
<img src="/images/3hff_pdb_example.png" width="300">
 
note: if you encounter issues visualizing your own pdb file, please contact us. Some pdb files need further clean-up before using directly in SCV.

#### Command line:
Rest API:
```shell script
curl --form 'psms=EQNEASPTPR,YCQEQDMCCR,ELAPGLHLR,GVVSDNCYPFSGR,C[143]TCHEGGHWECDQEPCLVDPDMIK,GRADECALPYLGATCYCDLFCN[115]R,GTNECDIETFVLGVWGR,EQNEASPTPR,GNYGWQAGN[115]HSAFWGMTLDEGIR,CPNGQVDSNDIYQVTPAYR,DLSWQVRSLLLDHNR,CNCALRPLCTWLR,RPGSRNRPGYGTGYF,RPDGDAASQPRTPILLLR,QSLRQELYVQDYASIDWPAQR,GTNGSQIWDTSFAIQALLEAGAHHR,ETLNQGLDFCRRKQR,SYFTDLPKAQTAHEGALN[115]GVTFYAK,CDGEANVFSDLHSLRQFTSR,ETFHGLKELAFSYLVWDSK,IKNIYVSDVLNMK' --form 'ptms={"N[115]":[0,255,8],"C[143]":[255,0,247]}' --form 'background_color=16777215' --form 'species=mouse' -X POST http://pepchem.org:35106/job
```
After sending a post request using curl, you will receive a job number which you can navigate to the site by going to http://pepchem.org.35106/view?job= returned job number from curl API

e.g. http://pepchem.org:35106/view?job=545cb66c-db8e-44f9-b495-454f25df5bec


### SCV local:
Please see example in [python/main.py](https://github.com/Gaolaboratory/SCV/tree/master/python/main.py)

1.Prepare your peptides in a list (with our without PTMs)

If you use Msfragger/Fragpipe as your database search tool, you could directly use the functions in [commons.py](http://pepchem.org:35091/blackjack/scv_local/blob/master/commons.py)
```python
peptide_list = commons.modified_peptide_from_psm('D:/data/native_protein_digestion/12072021/control/0240min/psm.tsv')
```

2.Prepare protein sequence dictionary (input peptides will be mapped on)

Directly from proteome fasta file:
```python
fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
protein_dict = commons.fasta_reader(fasta_file)
```

3.Define the PTM regex pattern and corresponding RGB color (you can skip this step if you don't want to show PTMs)
```python
ptm_color_dict = {'n\[43\]': [0, 0, 256]}  # here we show n terminal acetylation in blue
```
4.UniprotID for your protein of interest
```python
uniprot_id = 'P61956'
```

5.PDB file path for the protein above, if you don't have a PDB structure, you can predict it from [Alphafold2](https://github.com/deepmind/alphafold) or [RoseTTAFold](https://robetta.bakerlab.org/)
```python
pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61956-F1-model_v1.pdb'
```

6.HTML output path
```python
html_output = 'C:/Users/gao lab computer/PycharmProjects/SCV_local/P61956_test_ptm.html'
```

7.Call main function to generate 3d coverage in HTML
```python
_main_(peptide_list,
           protein_dict,
           uniprot_id,
           pdb_file=pdb_path,
           ptm_color_dict=ptm_color_dict,
           base_path=html_output)
```
### Result from html:

<img src="/images/example.png" width="250">

### How to cite:
- [Journal of Proteome Research](https://doi.org/10.1021/acs.jproteome.2c00358)
- Shao X, Grams C, Gao Y. Sequence Coverage Visualizer: A Web Application for Protein Sequence Coverage 3D Visualization. J Proteome Res. 2023 Feb 3;22(2):343-349. doi: 10.1021/acs.jproteome.2c00358. Epub 2022 Dec 13. PMID: 36511722; PMCID: PMC10232130.

For other tools developed by the Gao lab, see our website [https://lab.gy/](https://lab.gy/)

Or follow our twitter at https://twitter.com/gao_lab to see more exciting informatics tools!