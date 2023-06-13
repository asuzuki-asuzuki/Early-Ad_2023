# Haplotype assignment of mutations
## point mutations
1. Prepare a list of mutations (chr, position).
2. Extract information on the HP tags, PS tags, read names of mutations from Whatshap-happlotagged BAM files.  
`samtools mpileup -l mutation.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME happlotagged.bam > mutation.pileup`
3. Count reference and variant tags of each mutation position in each haplotype.


## SVs
