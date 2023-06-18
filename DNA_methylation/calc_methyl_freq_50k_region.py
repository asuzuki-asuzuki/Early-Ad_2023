import pandas as pd


def calc_50k_region_methylation(samplename:str) -> None:
    in_path = "~/methylation_frequency_{}.tsv".format(samplename)
    out_path = "~/50k_methyl_freq_{}_cut_off.tsv".format(samplename)
    chr_end_path = "~/hg38.chrom.sizes"

    df = pd.read_table(in_path)
    df_chr_end = pd.read_table(chr_end_path, header=None,index_col=0)

    sample_name = samplename
    pre_divided_site = df.iloc[0,1]//50000
    divided_block_count = 0
    methylated_site_number = 0
    read_number = 0
    methylated_read_number = 0
    with open(out_path, "w") as f:
        print("chromosome","start","end","sample_name","methylated_site_number","read_number","methylated_read_number","methylation_frequency",sep='\t',file=f)
        for i in df.itertuples():
            ref_chroms = []
            for j in range(1, 23):
                ref_chroms.append('chr' + str(j))
            ref_chroms.append('chrX')
            ref_chroms.append('chrY')
            if i[1] not in ref_chroms:
                continue
            
            if i[5] < 5:
                continue

            divided_site = i[2]//50000
            if i[2]%50000 == 0:
                divided_site = divided_site - 1

            if pre_divided_site == divided_site:
                methylated_site_number = methylated_site_number + i[4]
                read_number = read_number + i[5]
                methylated_read_number = methylated_read_number + i[6]
                chr_number = i[1]
                start_pos = divided_site * 50000 + 1
                end_pos = (divided_site + 1) * 50000
            else:
                methylated_frequency = methylated_read_number / read_number
                divided_block_count = divided_block_count + 1
                if i[1] != chr_number and df_chr_end.loc[chr_number,1] - start_pos < 50000:
                    end_pos = df_chr_end.loc[chr_number,1]
                print(chr_number,start_pos,end_pos,sample_name,methylated_site_number,read_number,methylated_read_number,methylated_frequency,sep='\t',file=f)
                pre_divided_site = divided_site
                methylated_site_number = 0 + i[4]
                read_number = 0 + i[5]
                methylated_read_number = 0 + i[6]
                chr_number = i[1]
                start_pos = divided_site * 50000 + 1
                end_pos = (divided_site + 1) * 50000
        methylated_frequency = methylated_read_number / read_number
        divided_block_count = divided_block_count + 1
        if df_chr_end.loc[chr_number,1] - start_pos < 50000:
            end_pos = df_chr_end.loc[chr_number,1]
        print(chr_number,start_pos,end_pos,sample_name,methylated_site_number,read_number,methylated_read_number,methylated_frequency,sep='\t',file=f)





def main():
    calc_50k_region_methylation(samplename)


if __name__ == "__main__":
    main()
