Fichiers source de l'annotation:
swissprot version 24 : <gitclonedir>/annotations/annot_source/uniprot_sprot.dat et <gitclonedir>/annotations/annot_source/uniprot_sprot.xml, uniprot_sprot.dat
metacyc version 20.1 : dossier 20.1 dans <gitclonedir>/annotations/annot_source/

Fichier de reference des EC numbers et des GO entries:
perl 8_write_reference_files.pl <gitclonedir> <version du dossier de metacyc> 8_intersection.out
output files:
<gitclonedir>/annotations/annot_files/SwissProt_ec.sql
<gitclonedir>/annotations/annot_files/SwissProt_go.sql
<gitclonedir>/annotations/annot_files/MetaCyc_ec.sql

<gitclonedir>/annotations/annot_files/SwissProt_MetaCyc_ec.sql
<gitclonedir>/annotations/annot_files/SwissProt_MetaCyc_go.sql

<gitclonedir>/annotations/annot_files/SwissProt_MetaCyc_ec.fa
<gitclonedir>/annotations/annot_files/SwissProt_MetaCyc_go.fa

Fichier de reference des fonctions swissprot (parsing uniprot_sprot.xml)
python 10A_parse_function_uniprot.py
output file:
<gitclonedir>/annotations/annot_files/swissprot.fa 

Fichier de reference des fonctions metacyc
perl 8_parse_proteins_dat.pl <gitclonedir> 20.1 > annot_files/parsing_proteins_metacyc
perl 8_parse_proteins_links_dat.pl <workdir> <metacyc version> > annot_files/metacyc_proteins_links> 
perl 10B_write_functions_references.pl <gitclonedir> 
output files:                    
annot_files/metacyc.fa 
                                          
                     
# group_annotation (EC numbers and Gene Ontology entries)
11A_run_group_hmm.bash
11B_hmm2ec.pl
11C_intersection_GO.pl
11D_list_group_function.pl 

GetSeq.pm                 
patch_none_metacyc_function.pl  

