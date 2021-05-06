## Phylogenetics

All phylogenies were constructed with IQTREE 2 (Minh et al. 2020) with 1000 UF Bootstraps (Hoang et al. 2018) and ModelFinder (Kalyaanamoorthy et al. 2017). 
```
iqtree2 -s [GENE].fasta --seqtype DNA --alrt 1000 -B 1000 -m MFP
```

Partitioned analyses (for genes and introns) were perfomed as follows: 
```
iqtree2 -s [GENE].fasta --seqtype DNA --alrt 1000 -B 1000 -m MFP -p partition_file.nex 
```

SNP phylogenies were constructed as above. The high-Fst data had no constant sites and we therefore ran ModelFinder to account for ascertainment bias as follows. 
```
iqtree2 -s high_fst_snps.fasta --seqtype DNA --alrt 1000 -B 1000 -m MFP+ASC 
```

<b> References: </b> 

IQTREE2: 

Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., in press. https://doi.org/10.1093/molbev/msaa015

ModelFinder: 

Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler, and Lars S Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. Nature Methods, 14:587â€“589. https://doi.org/10.1038/nmeth.4285

Ultra-Fast Bootstraps:  

Diep Thi Hoang, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh, and Le Sy Vinh (2018) UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol., 35: https://doi.org/10.1093/molbev/msx281
