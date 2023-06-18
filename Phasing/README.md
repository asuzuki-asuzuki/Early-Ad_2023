# Assignment of mutations to haplotypes
## point mutations (detected by GATK Mutect 2)
1. Prepare a list of mutations (chr, 0-based position; TSV format).  
[list.pl](./list.pl)  
2. Extract information on the HP tags, PS tags and read names of mutation positions from the Whatshap haplotagged BAM file.  
`samtools mpileup -l SNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME haplotagged.bam > SNV.pileup`  
`samtools mpileup -l indel.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME haplotagged.bam > indel.pileup`  
`samtools mpileup -l MNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME haplotagged.bam > MNV.pileup`  
4. Prepare a list of mutations with phased block information (chr, 0-based position, ref base, alt base, kind of mutations, phased block number).  
[list_phased.pl](./list_phased.pl)  
5. Count reference and variant tags of each mutation position in each haplotype.  
[extract_HP_SNV.pl](./extract_HP_SNV.pl), [extract_HP_indel.pl](./extract_HP_indel.pl), [extract_HP_MNV_qname.pl](./extract_HP_MNV_qname.pl), [extract_HP_MNV.pl](./extract_HP_MNV.pl)  
6. Determine the haplotype.  
[hantei_HP_SNV.pl](./hantei_HP_SNV.pl), [hantei_HP_indel.pl](./hantei_HP_indel.pl), [hantei_HP_MNV.pl](./hantei_HP_MNV.pl)  

## SVs (detected by nanomonsv)
1. Extract all read names assigned to the haplotypes (haplotagged.bam: BAM file with Haplotype (HP) tags determined by WhatsHap).
[extract_phased_readname.py](./extract_phased_readname.py)  
`python extract_phased_readname.py haplotagged.bam phased_read.txt`
3. Extract information on SVs (nanomonsv.result.filt.svtype.txt: output file of nanomonsv's post_filter.py/sv_type.py; nanomonsv.supporting_read.txt: output file of "nanomonsv get" command).
`python filter_supporting_read.py nanomonsv.result.filt.svtype.txt nanomonsv.supporting_read.txt nanomonsv.supporting_read_filt.txt`  
