# Assignment of mutations to haplotypes
## point mutations
1. Prepare a list of mutations (chr, 0-based position; TSV format).  
[list.pl](./list.pl)  
2. Extract information on the HP tags, PS tags, read names of mutations from Whatshap-happlotagged BAM files.  
`samtools mpileup -l SNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME happlotagged.bam > SNV.pileup`  
`samtools mpileup -l indel.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME happlotagged.bam > indel.pileup`  
`samtools mpileup -l MNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME happlotagged.bam > MNV.pileup`  
4. Prepare a list of mutations with phased block information (chr, 0-based position, ref base, alt base, kind of mutations, phased block number).  
[list_phased.pl](./list_phased.pl)  
5. Count reference and variant tags of each mutation position in each haplotype.  
[extract_HP_SNV.pl](./extract_HP_SNV.pl), [extract_HP_indel.pl](./extract_HP_indel.pl), [extract_HP_MNV_qname.pl](./extract_HP_MNV_qname.pl), [extract_HP_MNV.pl](./extract_HP_MNV.pl)  
6. Determining the haplotype.  
[hantei_HP_SNV.pl](./hantei_HP_SNV.pl), [hantei_HP_indel.pl](./hantei_HP_indel.pl), [hantei_HP_MNV.pl](./hantei_HP_MNV.pl)  

## SVs
