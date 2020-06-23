## Admixture Analysis 

We used the program Ohana to analyze admixture between <i> Plexaura homomalla </i> and <i> P. kukenthali </i>. 

1. Convert .vcf to .ped

<tt> plink --vcf p_homomalla_snps.vcf --recode 12 tab --geno 0.0 --out p_homomalla_snps --allow-extra-chr </tt> 

2. Convert .pdg to .dgm 

<tt> ped2dgm p_homomalla_snps.ped p_homomalla_snps.dgm </tt>

3. Generate q and f matrices

<tt> qpas p_homomalla_snps.dgm -k 2 -qo k2_q.matrix -fo k2_f.matrix -mi 100 </tt>

<tt> qpas p_homomalla_snps.dgm -k 3 -qo k3_q.matrix -fo k3_f.matrix -mi 100 </tt> 

4. Use the q-matrix to plot the likelihood for each individual of ancestry in each group.

<tt> python plot_q.py k2_q.matrix k2_bar_plot.pdf </tt>

<tt> python plot_q.py k3_q.matrix k3_bar_plot.pdf </tt>

<i> Note: After these pdfs were created, we remade the graphs in Adobe Illustrator. </i>

We ran this analysis for the whole SNP set and the high Fst SNPs. 

<b> Citations: </b>

Cheng, J.Y., T. Mailund, and R. Nielsen. 2017. Fast admixture analysis and population tree estimation for SNP and NGS data. <i> Bioinformatics </i> <b> 33(14): </b>2148-2155. 

Purcell, S., B. Neale, K. Todd-Brown, L. Thomas, M.A.R. Ferreira, D. Bender, J. Maller, P. Sklar, P.I.W. de Bakker, M.J. Daly, and P.C. Sham. 2007. PLINK: a toolset for whole-genome association and population-based linkage analysis. <i> American journal of Human Genetics </i> <b> 81(3): </b> 559-575. 
