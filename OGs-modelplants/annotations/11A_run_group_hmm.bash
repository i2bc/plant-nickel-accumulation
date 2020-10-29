for i in `seq n m`       
	do 
	muscle -in  <clonegitdir>/mario/TMP_METAPP_metOG_<rundate>/FINAL_GROUPS/$i.fa -out  <clonegitdir>/annotations/group_align/$i.fa ; 
	<clonegitdir>/mario/SCRIPTS/fasta2stockholm.pl  <clonegitdir>/annotations/group_align/$i.fa > <clonegitdir>/annotations/group_align/$i.sto ; 
	hmmbuild --cpu 25 <clonegitdir>/annotations/group_hmm/$i.msf <clonegitdir>/annotations/group_align/$i.sto ; 
	hmmsearch --cpu 25 -o /dev/null -E 1 --tblout <clonegitdir>/annotations/group_EC/EC$i.hmmsearch <clonegitdir>/annotations/group_hmm/$i.msf <clonegitdir>/annotations/annot_files/SwissProt_MetaCyc_ec.fa ; 
	hmmsearch --cpu 25 -o /dev/null -E 1 --tblout <clonegitdir>/annotations/group_GO/GO$i.hmmsearch  <clonegitdir>/annotations/group_hmm/$i.msf <clonegitdir>/annotations/annot_files/SwissProt_MetaCyc_go.fa ; 
	hmmsearch --cpu 25 -o /dev/null -E 1.e-6 --tblout <clonegitdir>/annotations/group_function/MC$i.hmmsearch <clonegitdir>/annotations/group_hmm/$i.msf <clonegitdir>/annotations/annot_files/metacyc.fa  
	hmmsearch --cpu 25 -o /dev/null -E 1.e-6 --tblout <clonegitdir>/annotations/group_function/SP$i.hmmsearch <clonegitdir>/annotations/group_hmm/$i.msf <clonegitdir>/annotations/annot_files/swissprot.fa
	done
