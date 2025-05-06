contacts="/home/anikolskaya/Mock_for_heat_1_arabi_grid.intersected.N2.for_bardic.bed"
annot="/home/anikolskaya/Arabidopsis_thaliana_genes_canonic_chrs_TAIR10_1_cleaned.for_bardic.bed"
chrsizes="/home/anikolskaya/GCA_000001735.2_TAIR10.1_genomic.fna.sizes"
pc="/home/anikolskaya/pc.txt"

bardic run $contacts $annot $chrsizes $pc ./results --qval_threshold 1 --cores 10 --top-percent 10