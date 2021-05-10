## Bioinformatic Pipeline for SNP Calling

Pipeline used to call Single Nucleotide Polymorphisms (SNPs) from <i> Plexaura homomalla </i> and <i> P. kukenthali </i> transcriptomes. Cleaned reads can be found in SRA under BioProject PRJNA675782, <i> de novo </i> assemblies can be found in this repository under "Data".



## 1. Download and assess <i> de novo </i> reference transcriptomes
 Use Benchmarking Universal Single-Copy Orthologs (BUSCO) to assess the completeness of each transcriptome. 

<tt> run_busco.py -i HF4_Trinity.assembly.fa -o HF4_Trinity_BUSCO -l eukaryota_odb9 -m tran </tt>

<tt> run_busco.py -i KG1_Trinity.assembly.fa -o KG1_Trinity_BUSCO -l eukaryota_odb9 -m tran </tt>

 Basic Usage:

<tt> run_busco.py -i [assembly_name.fa] -o [output_name] -l [lineage: either eukaryota_odb9 or metazoa_odb9] -m [mode: transcriptome] </tt>

<b> Citation: </b> Simao, F.A., R.M. Waterhouse, P. Ioannidis, E.V. Kriventseva, and E.M. Zdobnov. 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. <i> Bioinformatics </i> <b> 31(19) </b> 3210-3212.  

## 2. Map reads to reference transcriptome 
Use BowTie2 to map reads from each colony to the reference transcriptome (HF4, based on higher completeness). 

<b> First, </b> build an index using BowTie-build- this will be the refrerence to which reads are mapped. 

<tt> bowtie2-build -f HF4_Trinity.assembly.fa HF4_Trinity.assembly </tt>

 Basic usage:

<tt> bowtie2-build -f [input_assembly.fa] [index_name] </tt> 

<b> Second, </b> use BowTie2 to map all reads to the reference transcriptome. The output alignments will be in the .sam format. In this case, we did not use paired reads and instead only use the forward (_1) reads.  

<tt> bowtie2 -x HF4_Trinity.assembly -U HF1_R1_paired_clean.fq -S HF1.sam --no-unal --threads 6; </tt>

<tt> bowtie2 -x HF4_Trinity.assembly -U HF2_R1_paired_clean.fq -S HF2.sam --no-unal --threads 6; </tt>

<tt> bowtie2 -x HF4_Trinity.assembly -U HF3_R1_paired_clean.fq -S HF3.sam --no-unal --threads 6; </tt> 

<tt> bowtie2 -x HF4_Trinity.assembly -U HF4_R1_paired_clean.fq -S HF4.sam --no-unal --threads 6; </tt> 

<tt> bowtie2 -x HF4_Trinity.assembly -U HF5_R1_paired_clean.fq -S HF5.sam --no-unal --threads 6; </tt>

<tt> bowtie2 -x HF4_Trinity.assembly -U KG1_R1_paired_clean.fq -S KG1.sam --no-unal --threads 6; </tt> 

<tt> bowtie2 -x HF4_Trinity.assembly -U KG3_R1_paired_clean.fq -S KG3.sam --no-unal --threads 6; </tt>

<tt> bowtie2 -x HF4_Trinity.assembly -U KG4_R1_paired_clean.fq -S KG4.sam --no-unal --threads 6 </tt> 

Basic usage: 

<tt> bowtie2 -x [index_name] -U [individual_1_unpaired_reads] -S [individual_1_SAM_file] --no-unal --threads [# threads] </tt> 

<b> Citations: </b>

Langmead, B. and S.L. Salzburg. 2012. Fast gapped-read alignment with Bowtie2. <i> Nature Methods </i> <b> 9(4): </b> 357-359. 

Langmead, B., C. Wilks, V. Antonescu, and R. Charles. 2018. Scaling read aligners to hundreds of threads on general-purpose processors. <i> Bioinformatics </i> <b> 35(3): </b> 421-432. 

## 3. Convert .sam to .bam
Use Samtools to convert .sam alignment files into .bam files. The output will be a binary version of the .sam alignments. 

<tt> samtools view -S -b HF1.sam > HF1.bam -@ 6; samtools view -S -b HF2.sam > HF2.bam -@ 6; </tt> 

<tt> samtools view -S -b HF3.sam > HF3.bam -@ 6; samtools view -S -b HF4.sam > HF4.bam -@ 6; </tt> 

<tt> samtools view -S -b HF5.sam > HF5.bam -@ 6; samtools view -S -b KG1.sam > KG1.bam -@ 6; </tt> 

<tt> samtools view -S -b KG3.sam > KG3.bam -@ 6; samtools view -S -b KG4.sam > KG4.bam -@ 6; </tt> 

Basic Usage:

<tt> samtools view -S -b [input_sam_file] > [output_bam_file] -@ [# threads] </tt> 

<b> Citation: </b> Li, H., B. Handsaker, A. Wysoker, T. Fenell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R. Durbin, and 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map forat and SAMtools. <i> Bioinformatics </i> <b> 25(16): </b> 2078-2079. 

## 4. Sort and index .bam files 
Use Samtools to sort and then index resulting .bam files.

<b> First, </b> sort the .bam files.

<tt> samtools sort -@6 HF1.bam -o HF1_sorted.bam; samtools sort -@6 HF2.bam -o HF2_sorted.bam </tt>

<tt> samtools sort -@6 HF3.bam -o HF3_sorted.bam; samtools sort -@6 HF4.bam -o HF4_sorted.bam </tt>

<tt> samtools sort -@6 HF5.bam -o HF5_sorted.bam; samtools sort -@6 KG1.bam -o KG1_sorted.bam </tt>

<tt> samtools sort -@6 KG3.bam -o KG3_sorted.bam; samtools sort -@6 KG4.bam -o KG4_sorted.bam </tt>

Basic Usage: 

<tt> samtools sort -@ [# threads] [input_bam_file] -o [sorted_bam_file] </tt>

<b> Second, </b> index the sorted bam files.

<tt> samtools index HF1_sorted.bam -@6; samtools index HF2_sorted.bam -@6 </tt> 

<tt> samtools index HF3_sorted.bam -@6; samtools index HF4_sorted.bam -@6 </tt> 

<tt> samtools index HF5_sorted.bam -@6; samtools index KG1_sorted.bam -@6 </tt> 

<tt> samtools index KG3_sorted.bam -@6; samtools index KG4_sorted.bam -@6 </tt> 

Basic Usage:

<tt> samtools index [sorted_bam_file] -@ [# threads] </tt> 

<b> Citation: </b> Li, H., B. Handsaker, A. Wysoker, T. Fenell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R. Durbin, and 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map forat and SAMtools. <i> Bioinformatics </i> <b> 25(16): </b> 2078-2079. 

## 5. Variant Calling
Use bcftools to call variants and SNPs. Note that many of these parameters should be adjusted for individual cases. 

<b> First, </b> use the <tt> mpileup </tt> command. 

<tt> bcftools mpileup -f HF4_Trinity.assembly.fa -b sample_ids.txt -o mpileup_out.vcf -q 5 -a AD, GT, PL, -Ov --threads 7 -I -C 50 -Q25 </tt>

Basic Usage:

<tt> bcftools mpileup -f [indexed_ref_fasta] -b [file_names.txt] -o [output_file] -q [mapping quality, 1 to 5] -a [annotations] -Ov [output is VCF] --threads [# threads] -I [exclude indels] -C [rec. value is 50]  -Q [min. base quality] </tt>

<b> Second, </b> use the <tt> call </tt> command.

<tt> bcftools call mpileup_out.vcf -Ov -o call_out.vcf --threads 7 -V indels -v -c </tt>

Basic Usage:

<tt> bcftools call [input_vcf] -Ov [output is VCF] -o [output_file.vcf] --threads [# threads] -V [indels, exclude indels] -v [only output variant sites] -c [use consensus caller] </tt>

<b> Citation: </b> Li, H. 2011. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. <i> Bioinformatics </i> <b> 27(1): </b> 2987-2993. 

## 6. Filter SNPs

<b> First, </b> Use Samtools and VCFTools to filter the vcf file from #5. 

<b> First, </b> use the perl script <tt> vcfutils.pl </tt> in samtools to filter for SNP depth. Note that confident SNP calling in heterozygous genomes requires at least 13-15X coverage (Meynert et al. 2013, Song et al. 2016).    

<tt> vcfutils.pl VarFilter -d 30 call_out.vcf > vcf_utils_filt.vcf </tt> 

Basic Usage: 

<tt> vcfutils.pl VarFilter -d [integer, x coverage] [input.vcf] > [output.vcf] </tt>

<b> Second, </b> use VCFTools for further filtering. 

<tt> vcftools --vcf vcf_utils_filt.vcf --max-missing 1 --thin 5000 --recode </tt> 

Basic Usage:

<tt> vcftools --vcf [input.vcf] --max-missing [0-1, 1= no missing data allowed] --thin [integer, only keep SNPs xbp apart] --recode [writes new .vcf] </tt> 

<b> Citation: </b> Danecek, P., A. Auton, G. Abecasis, C.A. Albers, E. Banks, M.A DePristo, R. Handsaker, G. Lunter, G. Marth, S.T. Sherry, G. McVean, R. Durbin, and 1000 Genomes Project Analysis Group. 2011. The variant call format and VCFTools. <i> Bioinformatics </i> <b> 27(15): </b> 2156-2158. 
