# mt-paths
Hacky program to find a set of best "paths" through an mtDNA phylogeny, given a set of aligned reads.

## current example usage:

```
python bin/mt-paths.py data/mtDNA_tree_Build_17.csv -m 6776 9055 10550 12372 14167 14766 -r /mnt/expressions/jonas_platzek/2_references/Homo_sapiens.fasta -b /mnt/expressions/jonas_platzek/0_playground/RGK/bam/Cap.E.8300.Hominidae.Homo_sapiens_deduped.bam --genbank /mnt/expressions/jonas_platzek/0_playground/RGK/genbank_files/NC_012920.1.gb
```
