We downloaded and store the fasta formatted proteomes in <path_to_OGs-modelplants>/modelsfasta/ 
 
# blast step
cd <path_to_OGs-modelplants>/run_blast
perl allblast_with_threads.pl <path_to_OGs-modelplants>/modelsfasta/ <path_to_OGs-modelplants>/modelsblast/ <cpu>
output files in <path_to_OGs-modelplants>/modelsblast/

# numbering step
cd <path_to_OGs-modelplants>/manage_data
perl 1_numerotation_genomes.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/modelsfasta/ 
perl 2_numerotation_blast.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/modelsblast/ 
output files in <path_to_OGs-modelplants>/numfast and <path_to_OGs-modelplants>/numblast



# BRH run
cd <path_to_OGs-modelplants>/BRH
perl 1_BRH.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numblast/ <cpu> 
perl 2_BRH.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ <cpu> 
perl 3_BRH.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ 
output files: STEP3/brh_orthologous_groups.tab

# Inparanoid run
cd <path_to_OGs-modelplants>/Inparanoid
perl 1_parser_pour_Inparanoid.pl 40 <path_to_OGs-modelplants>/numblast/ <cpu> 
perl 2_inparanoid.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ <cpu> 
perl 3_inparanoid.pl <path_to_OGs-modelplants>/ > inparanoid.tab
output files: STEP3/inparanoid_inparalogs.sql and STEP3/inparanoid_orthologous_groups.sql

# Phylogeny run 
cd <path_to_OGs-modelplants>/Phylogeny
perl 1_phylogeny_launch_jobs.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ <path_to_OGs-modelplants>/numblast/ <cpu> 
perl 2_phylogeny.pl 70 <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ 
perl 3_phylogeny_launch_jobs.pl 12 <cpu> <path_to_OGs-modelplants>/ 
perl 4_phylogeny_launch_jobs.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ <cpu> 
perl 5_phylogeny_modif.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/numfasta/ 
output files: STEP5/phylogeny_inparalogs.sql  and STEP5/phylogeny_orthologous_groups.sql

# OrthoFinder run (use precalculated blasts via symbolic links to fasta and blast directories)
cd marioplants/OrthoFinder/
python2 orthofinder.py -b input
input files: 
SpeciesIDs.txt
SequenceIDs.txt
Symbolic links to fasta and precalculated blasts directories
output files:
clusters_OrthoFinder_v0.4.0_I1.5.txt
clusters_OrthoFinder_v0.4.0_I1.5.txt_id_pairs.txt
OrthoFinder_v0.4.0_graph.txt
OrthologousGroups.csv
OrthologousGroups.txt
OrthologousGroups_UnassignedGenes.csv

# formatting output files for MaRiO
cd <path_to_OGs-modelplants>/scripts
perl phylogeny2mario.pl <path_to_OGs-modelplants>/
perl inparanoid2mario.pl <path_to_OGs-modelplants>/
perl brh2mario.pl <path_to_OGs-modelplants>/
perl orthofinder2mario.pl <path_to_OGs-modelplants>/
output files: mario_converted/brh_groups, mario_converted/inparanoid_groups, mario_converted/orthofinder_groups, mario_converted/phylogeny_groups
perl back_to_name-list.pl ../mario_input/brh_groups ';' ../OrthoFinder/input/SpeciesIDs.txt': ' ../numfiles/ ': ' 
perl back_to_name-list.pl ../mario_input/inparanoid_mario_groups ';' ../OrthoFinder/input/SpeciesIDs.txt ': ' ../numfiles/ ': '
output files: mario_input/brh_groups, mario_input/inparanoid_groups, mario_input/orthofinder_groups, mario_input/phylogeny_groups

# mario
cd <path_to_OGs-modelplants>/mario
perl mario.pl -i <path_to_OGs-modelplants>/mario_input/ -f <path_to_OGs-modelplants>/modelsfasta/ -cpu <cpu>
output (fasta formatted files): <path_to_OGs-modelplants>/mario/RESULTS_metOG_7_10_2017/FINAL_GROUPS.bz2 

We download and store in the <path_to_OGs-modelplants>/annotations/annot_source/ the thereafter ressources:
- from Swiss-Prot version 24 (ftp:/ftp.uniprot.org) : <path_to_OGs-modelplants>/annotations/annot_source/uniprot_sprot.dat et <path_to_OGs-modelplants>/annotations/annot_source/uniprot_sprot.xml
- from MetaCyc version 20.1 (https:/metacyc.org) : <path_to_OGs-modelplants>/annotations/annot_source/20.1

# group_annotation (function, EC numbers and Gene Ontology entries)

## Building reference files for EC numbers and GO entries:
perl 8_write_reference_files.pl <path_to_OGs-modelplants>/ <version du dossier de metacyc> 8_intersection.out
output files:
<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_ec.sql
<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_go.sql
<path_to_OGs-modelplants>/annotations/annot_files/MetaCyc_ec.sql

<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_MetaCyc_ec.sql
<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_MetaCyc_go.sql

<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_MetaCyc_ec.fa
<path_to_OGs-modelplants>/annotations/annot_files/SwissProt_MetaCyc_go.fa

## Building reference files for swissprot functions (parsing uniprot_sprot.xml)
python 10A_parse_function_uniprot.py
output file:
<path_to_OGs-modelplants>/annotations/annot_files/swissprot.fa 

## Building reference files for metacyc functions:
cd <path_to_OGs-modelplants>/annotations
perl 8_parse_proteins_dat.pl <path_to_OGs-modelplants>/ 20.1 > annot_files/parsing_proteins_metacyc
perl 8_parse_proteins_links_dat.pl <workdir> <metacyc version> > annot_files/metacyc_proteins_links> 
perl 10B_write_functions_references.pl <path_to_OGs-modelplants>/ 
output files:                    
<path_to_OGs-modelplants>/annotations/annot_files/metacyc.fa 

perl 9_intersect_ref_fungipath.pl <path_to_OGs-modelplants>/ ../liste_plants <path_to_OGs-modelplants>/modelsfasta/ 
bash 11A_run_group_hmm.bash
perl 11B_hmm2ec.pl <path_to_OGs-modelplants>/
perl 11C_intersection_GO.pl <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/ <path_to_OGs-modelplants>/mario/TMP_METAPP_metOG_7_10_2017/FINAL_GROUPS/
perl 11D_list_group_function.pl <path_to_OGs-modelplants>/
perl fill_groupes_table.pl > bestfunction.tab









