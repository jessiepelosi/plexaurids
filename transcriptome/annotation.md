## Annotating Octocoral Transcriptomes 

1. Make blast databases from peptides downloaded for Cnidarians in Ensembl Metazoa 47, remaining taxa in Ensembl Metazoa 47, and Cnidarians in Uniprot (SwissProt and Trembl). 

<tt> makeblastdb -in cnidarians_ensembl.fasta -parse_seqikds -dbtype prot -o cnidarian_ensembl_db </tt>

<tt> makeblastdb -in metazoa_ensembl.fasta -parse_seqikds -dbtype prot -o metazoa_blast_db </tt>

<tt>  makeblastdb -in uniprot_cnidaria.fasta -parse_seqikds -dbtype prot -o uniprot_cnidaria_swiss_trembl_db </tt>


2. Annotate transcriptomes using <i> blastx </i> with Cnidarians in Ensebml Metazoa 47. 

<tt> blastx -db cnidarian_ensembl_db -query HF4_Trinity.fasta -out HF4_ensembl_cnidaria_annotation -evalue 0.00001 -max_target_seqs 10 </tt>

<tt> blastx -db cnidarian_ensembl_db -query KG1_Trinity.fasta -out KG1_ensembl_cnidaria_annotation -evalue 0.00001 -max_target_seqs 10 </tt>

3. Annotate transcriptomes using <i> blastx </i> with remaining taxa in Ensembl Metazoa 47. 

<tt> blastx -db metazoa_blast_db -query HF4_Trinity.fasta -out HF4_ensembl_metazoa_annotation -evalue 0.00001 -max_target_seqs 10 </tt>

<tt> blastx -db metazoa_blast_db -query KG1_Trinity.fasta -out KG1_ensembl_metazoa_annotation -evalue 0.00001 -max_target_seqs 10 </tt>

4. Annotate transcriptomes using <i> blastx </i> with Cnidarians from Uniprot. 

<tt> blastx -db uniprot_cnidaria_swiss_trembl_db -query HF4_Trinity.fasta -out HF4_uniprot_annotation -evalue 0.00001 -max_target_seqs 10 </tt> 

<tt> blastx -db uniprot_cnidaria_swiss_trembl_db -query KG1_Trinity.fasta -out KG1_uniprot_annotation -evalue 0.00001 -max_target_seqs 10 </tt> 

5. Use these annotations in the R script "annotation_filtering.R" to finalize the annotations. 

<b> Citations: </b> 

Altschul, S.F., W. Gish, W. Miller, E.W. Myers, D.J. Lipman. 1990. Basical local alignment search tools. <i> Journal of Molecular Biology </i> <b> 215: </b> 403-410. 

Cunningham, F., P. Achuthan, W. Akanni, et al. 2019. Ensembl 2019. <i> Nucleic Acids Research </i> <b> 47(D1): </b> D745-D751. 

The Uniprot Consortium. 2019. UniProt: A worldwide hub of protein knowledge. <i> Nucleic Acids Research </i> <b> 47: </b> D506-D515. 
