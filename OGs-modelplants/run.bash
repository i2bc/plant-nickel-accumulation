
We download and store in the OGs-modelplants/fasta directory the fasta formatted proteomes
We download and store in the OGs-modelplants/annotations/annot_source/ the thereafter ressources
swissprot version 24 (ftp://ftp.uniprot.org) : OGs-modelplants/annotations/annot_source/uniprot_sprot.dat et OGs-modelplants/annotations/annot_source/uniprot_sprot.xml, uniprot_sprot.dat
metacyc version 20.1 (https://metacyc.org) : OGs-modelplants/annotations/annot_source/20.1

Building reference files for EC numbers and GO entries:
perl 8_write_reference_files.pl OGs-modelplants <version du dossier de metacyc> 8_intersection.out
output files:
OGs-modelplants/annotations/annot_files/SwissProt_ec.sql
OGs-modelplants/annotations/annot_files/SwissProt_go.sql
OGs-modelplants/annotations/annot_files/MetaCyc_ec.sql

OGs-modelplants/annotations/annot_files/SwissProt_MetaCyc_ec.sql
OGs-modelplants/annotations/annot_files/SwissProt_MetaCyc_go.sql

OGs-modelplants/annotations/annot_files/SwissProt_MetaCyc_ec.fa
OGs-modelplants/annotations/annot_files/SwissProt_MetaCyc_go.fa

Building reference files for swissprot functions (parsing uniprot_sprot.xml)
python 10A_parse_function_uniprot.py
output file:
OGs-modelplants/annotations/annot_files/swissprot.fa 

Building reference files for metacyc functions:
perl 8_parse_proteins_dat.pl OGs-modelplants 20.1 > annot_files/parsing_proteins_metacyc
perl 8_parse_proteins_links_dat.pl <workdir> <metacyc version> > annot_files/metacyc_proteins_links> 
perl 10B_write_functions_references.pl OGs-modelplants 
output files:                    
annot_files/metacyc.fa 
                                          
                     
# group_annotation (EC numbers and Gene Ontology entries)
11A_run_group_hmm.bash
11B_hmm2ec.pl
11C_intersection_GO.pl
11D_list_group_function.pl 

GetSeq.pm                 
patch_none_metacyc_function.pl  

