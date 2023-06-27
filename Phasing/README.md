# Assignment of mutations to haplotypes
## point mutations (detected by GATK Mutect 2)
1. Prepare a list of mutations (chr, 0-based position; TSV format).  
[list.pl](./list.pl)  
2. Extract information on the HP tags, PS tags and read names of mutation positions from the Whatshap haplotagged BAM file using SAMtools (version 1.12).  
```
samtools mpileup  
      -l SNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME
      haplotagged.bam > SNV.pileup
```  
```
samtools mpileup  
      -l indel.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME  
      haplotagged.bam > indel.pileup
```  
```
samtools mpileup  
      -l MNV.list --min-BQ 5 --min-MQ 20 --excl-flags 2304 --output-extra HP,PS,QNAME  
      haplotagged.bam > MNV.pileup
```  
3. Prepare a list of mutations with phased block information (chr, 0-based position, ref base, alt base, kind of mutations, phased block number).  
[list_phased.pl](./list_phased.pl)  
4. Count reference and variant tags of each mutation position in each haplotype.  
[extract_HP_SNV.pl](./extract_HP_SNV.pl), [extract_HP_indel.pl](./extract_HP_indel.pl), [extract_HP_MNV_qname.pl](./extract_HP_MNV_qname.pl), [extract_HP_MNV.pl](./extract_HP_MNV.pl)  
5. Determine the haplotype.  
[hantei_HP_SNV.pl](./hantei_HP_SNV.pl), [hantei_HP_indel.pl](./hantei_HP_indel.pl), [hantei_HP_MNV.pl](./hantei_HP_MNV.pl)  

## SVs (detected by nanomonsv)
1. Extract all read names assigned to the haplotypes (**haplotagged.bam**: BAM file with Haplotype (HP) tags determined by WhatsHap).
[extract_phased_readname.py](./extract_phased_readname.py)  
```
python extract_phased_readname.py haplotagged.bam phased_read.txt
```  
2. Extract information on SVs (**nanomonsv.result.filt.svtype.txt**: output file of nanomonsv's post_filter.py/sv_type.py; **nanomonsv.supporting_read.txt**: output file of "nanomonsv get" command).  
[filter_supporting_read.py](./filter_supporting_read.py)  
```
python filter_supporting_read.py
      nanomonsv.result.filt.svtype.txt nanomonsv.supporting_read.txt nanomonsv.supporting_read_filt.txt
```
3. Add phasing information (**block-list.tsv**: WhatsHap output file for phased block regions).  
[add_phaseinfo_nanomonsv.py](./add_phaseinfo_nanomonsv.py)  
```
python add_phaseinfo_nanomonsv.py
      nanomonsv.result.filt.svtype.txt nanomonsv.supporting_read_filt.txt phased_read.txt  
      block-list.tsv nanomonsv.supporting_read_filt_phased.txt
```  
4. Divide reads into primary or supplementary alignment (**q20.bam**: BAM files filtered by MAPQ 20 or higher in this project).  
[extract_reads_nanomonsv.py](./extract_reads_nanomonsv.py)  
```
python extract_reads_nanomonsv.py  
      q20.bam nanomonsv.supporting_read_filt_phased.txt  
      nanomonsv_reads_primary.bam supplementary_tmp.bam
```  
5. Extract uniquely mapped reads from Supplementary alignment BAM file.  
[extract_uniq_reads.py](./extract_uniq_reads.py)  
```
python extract_uniq_reads.py  
      supplementary_tmp.bam nanomonsv_reads_supplementary.bam
```  
6. Sort and index BAM files.  
```
samtools sort -@ 10 nanomonsv_reads_primary.bam -o nanomonsv_reads_primary.sort.bam  
samtools sort -@ 10 nanomonsv_reads_supplementary.bam -o nanomonsv_reads_supplementary.sort.bam  
samtools index nanomonsv_reads_primary.sort.bam  
samtools index nanomonsv_reads_supplementary.sort.bam
```  
7. Extract SNP information from supplementary alignments (**whatshap.phased.vcf**: VCF file of germline SNPs with phase information from WhatsHap phase).  
```
cat  whatshap.phased.vcf  | grep -v "#" | cut -f 1,2,4,5,9,10 | grep PS > germlinesnps_pos_list.txt  
samtools mpileup -l germlinesnps_pos_list.txt --min-BQ 5 --min-MQ 20 --output-QNAME  
      nanomonsv_reads_supplementary.sort.bam > samtools_mpileup_list_supplementary.txt
```  
8. Assign supplementary alignment reads to haplotypes when ≥2 phased SNPs were detected and ≥70% of these SNPs were located on one haplotype.  
[phasing_supplementary_read.py](./phasing_supplementary_read.py)  
```
python phasing_supplementary_read.py  
      nanomonsv.supporting_read_filt.txt germlinesnps_pos_list.txt samtools_mpileup_list_supplementary.txt  
      block-list.tsv nanomonsv.result.filt.svtype.txt nanomonsv.supporting_read_supplementary_phased.txt
```  
9. Count phased reads for each SV breakpoint.  
[stat_phasing_sv.py](./stat_phasing_sv.py)  
```
python stat_phasing_sv.py  
      nanomonsv.supporting_read_filt_phased.txt nanomonsv.supporting_read_supplementary_phased.txt  
      nanomonsv_phasing_summary.txt
```  
10. Assign each SV to the haplotype when three or more supporting reads were assigned to the haplotype, and ≥70% of these supporting reads harbored one haplotype tag (HP1 or HP2).  
[stat_phasing_sv2_fix.py](./stat_phasing_sv2_fix.py)  
```
python stat_phasing_sv2_fix.py  
      nanomonsv_phasing_summary.txt nanomonsv_phasing_summary2_fix.txt block-list.tsv
```  

# Determination of the order of mutation occurrence on long reads  
1. Count base patterns of two closely located SNVs at a haplotype level (**SNV.PASS.vcf**: output file of GATK Mutect 2 and Annovar with filter PASS; **SNV.pileup**: output file of SAMtools mpileup for SNV position).  
[call_phased_snv.py](./call_phased_snv.py)  
```
python call_phased_snv.py SNV.PASS.vcf SNV.pileup > call_phased.tsv
```  
2. Extract mutation pairs of which the occurrence order can be inferred (threshold: set to 2 in this project).  
[count_ordered_snv.py](./count_ordered_snv.py)  
```
python count_ordered_snv.py call_phased.tsv threshold > position.tsv
```

