# COVID19Consensus
This plugin builds a consensus sequence of the SARS-CoV-2 virus present in the sample through two rounds of mapping and variant calling (only detecting mutations with allele frequency >= 51%), and masks low and no coverage positions.

## 2 rounds of align-call-refine cycle:
Reads are aligned to hg19_human_trxo_and_cor.fasta using TMAP. Variant positions in the 2019-nCoV are called using TVC with chip-specific ampliseq germline low stringency parameters except that variant allele frequency thresholds are changed to 0.51. Then 2019-nCoV in the reference is refined to represent the major allele at each variant site. This align-call-refine cycle is iterated twice, to minimize reference bias in the 2019-nCoV.

Additionally, it will also mask low (< user specified minimum coverage in the plugin interface) and no coverage positions. When calculating depth, it will take the mapping quality which can also be set in the plugin interface into account. 

Note Ion Ampliseq Coronavirus panel just covers the genomic region: 42-29842 (1-29903), so it will delete extra bases on the two ends ([1,41], [length-60,length]) in the final contig.

