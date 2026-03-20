################################################################################
# SETUP ENVIRONMENT

# metaxcan environment
sudo apt-get -y install git
git clone https://github.com/hakyimlab/MetaXcan
git clone https://github.com/hakyimlab/summary-gwas-imputation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sudo bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3
export PATH=$PATH:/opt/miniconda3/bin
conda config --remove channels defaults
conda config --add channels conda-forge
conda init
bash
conda activate imlabtools

# conda environment for harmonisation scripts
cd MetaXcan
git show 0b7c10d633d3d7bfe794e4f35bd8190356bb2514:software/conda_env.yaml > software/conda_env_old.yaml 
cd ~
conda env create -f MetaXcan/software/conda_env_old.yaml --yes

# "tmux" "tmux a 0" "tmux det"
sudo apt-get -y install tmux
echo "set -g mouse on" >> ~/.tmux.conf

# download required GWAS summary data
pip install gdown
gdown [[GWAS DATA]] --fuzzy
mkdir ~/data
tar -xvzf data.tar.gz -C data
rm data.tar.gz

sudo apt-get -y install unzip
sudo apt-get -y install bcftools

# variables
ROOT="[[PLACEHOLDER]]"
GWAS_TOOLS="$ROOT/summary-gwas-imputation/src"
METAXCAN="$ROOT/MetaXcan/software"
OUTPUT="$ROOT/outputs"
DATA="$ROOT/data"
MODEL="$DATA/models/eqtl/mashr"
LIFTOVER="$DATA/liftover"
mkdir $OUTPUT

################################################################################
# BUILD REFERENCE DATASETS
# make samples.txt from "GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS" based on RACE=3 and ETHNICTY=0
# using https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000424/phs000424.v8.p2/pheno_variable_summaries/phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml

# make gtex_v8_eur_filtered.txt.gz
echo -ne "varID\t" | gzip > $DATA/gtex_reference/gtex_v8_eur_filtered.txt.gz
bcftools view $DATA/gtex_reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf -S $DATA/gtex_reference/samples.txt --force-samples -Ou |  bcftools query -l | tr '\n' '\t' | sed 's/\t$/\n/' | gzip >> $DATA/gtex_reference/gtex_v8_eur_filtered.txt.gz

NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Starting at $NOW"
bcftools view -S $DATA/gtex_reference/samples.txt --force-samples $DATA/gtex_reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf -Ou | \
bcftools query -f '%ID[\t%GT]\n' | \
awk '
{
for (i = 1; i <= NF; i++) {
    if (substr($i,0,1) == "c") {
        printf("%s",$i)
    } else if ( substr($i, 0, 1) == ".") {
        printf("\tNA")
    } else if ($i ~ "[0-9]|[0-9]") {
        n = split($i, array, "|")
        printf("\t%d",array[1]+array[2])
    } else {
        #printf("\t%s",$i)
        printf("Unexpected: %s",$i)
        exit 1
    }
}
printf("\n")
}
' | gzip >> $DATA/gtex_reference/gtex_v8_eur_filtered.txt.gz
NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Ending at $NOW"

# make gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz
python3 $GWAS_TOOLS/get_reference_metadata.py \
	-genotype $DATA/gtex_reference/gtex_v8_eur_filtered.txt.gz \
	-annotation $DATA/gtex_reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz \
	-rsid_column rs_id_dbSNP151_GRCh38p7 \
	-filter MAF 0.01 \
	-filter TOP_CHR_POS_BY_FREQ \
	-output $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz

# make PARQUETs of gtex_v8_eur_filtered_maf0.01_monoallelic_variants
python3 $GWAS_TOOLS/model_training_genotype_to_parquet.py \
	-input_genotype_file $DATA/gtex_reference/gtex_v8_eur_filtered.txt.gz \
	-snp_annotation_file $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-parsimony 9 \
	--impute_to_mean \
	--split_by_chromosome \
	--only_in_key \
	-output_prefix genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants

################################################################################
# DATA HARMONISATION
# tmux
# conda activate imlabtools
# test headings order:  zcat filename.txt.gz | head -n 5

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $ROOT/pineiro_fjell_autosome_plink_snpsonly_rmdup_Didac.PHENO2.glm.linear.withNEA.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/pineiro_fjell_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map CHROM chromosome \
	-output_column_map POS position \
	-output_column_map ID variant_id \
	-output_column_map A1 effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	-output_column_map non_effect_allele non_effect_allele \
	--insert_value sample_size 35712 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $ROOT/kim_lee_2.txt \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/kim_lee_gwas_harmonised.txt.gz \
	--chromosome_format \
	-separator " " \
	-output_column_map CHR chromosome \
	-output_column_map POS position \
	-output_column_map SNPID variant_id \
	-output_column_map Allele2 effect_allele \
	-output_column_map AF_Allele2 frequency \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map p.value pvalue \
	-output_column_map Allele1 non_effect_allele \
	--insert_value sample_size 34129 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $ROOT/GCST90558034.tsv \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/zhao_zhao_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map chromosome chromosome \
	-output_column_map base_pair_location position \
	-output_column_map effect_allele effect_allele \
	-output_column_map beta effect_size \
	-output_column_map standard_error standard_error \
	-output_column_map p_value pvalue \
	-output_column_map other_allele non_effect_allele \
	--insert_value sample_size 30400 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $ROOT/kaufman_westlye_Brainage_GWAS_sumstat_final \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/kaufmann_westlye_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map CHR chromosome \
	-output_column_map SNP variant_id \
	-output_column_map BP position \
	-output_column_map A1 effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map PVAL pvalue \
	-output_column_map A2 non_effect_allele \
	--insert_value sample_size 20170 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/fan_huang/bag.txt.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/fan_huang_gwas_harmonised.txt.gz \
	--chromosome_format \
	-separator " " \
	-output_column_map CHR chromosome \
	-output_column_map SNP variant_id \
	-output_column_map BP position \
	-output_column_map A1 effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	-output_column_map A2 non_effect_allele \
	-output_column_map MAF frequency \
	--insert_value sample_size 31520 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/leonardsen_wang/BA_meta_cleaned.txt.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/leonardsen_wang_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map CHR chromosome \
	-output_column_map BP position \
	-output_column_map SNP variant_id \
	-output_column_map A1 effect_allele \
	-output_column_map A2 non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	-output_column_map FRQ frequency \
	--insert_value sample_size 28104 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/wen_davatzikos/bag_pheno_normalized_residualized.wm.glm.linear \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/wen_davatzikos_wm_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map "#CHROM" chromosome \
	-output_column_map POS position \
	-output_column_map ID variant_id \
	-output_column_map ALT effect_allele \
	-output_column_map REF non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 31129 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/wen_davatzikos/bag_pheno_normalized_residualized.gm.glm.linear \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/wen_davatzikos_gm_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map "#CHROM" chromosome \
	-output_column_map POS position \
	-output_column_map ID variant_id \
	-output_column_map ALT effect_allele \
	-output_column_map REF non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 31129 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/wen_davatzikos/bag_pheno_normalized_residualized.fc.glm.linear \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/wen_davatzikos_fc_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map "#CHROM" chromosome \
	-output_column_map POS position \
	-output_column_map ID variant_id \
	-output_column_map ALT effect_allele \
	-output_column_map REF non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 31129 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/wen_davatzikos/Brain_age_gap_PLINK_EUR.glm.linear \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/wen_davatzikos_combined_gwas_harmonised.txt.gz \
	--chromosome_format \
	-output_column_map "#CHROM" chromosome \
	-output_column_map POS position \
	-output_column_map ID variant_id \
	-output_column_map ALT effect_allele \
	-output_column_map REF non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 29400 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pminuslog10 pvalue zscore effect_size standard_error sample_size n_cases frequency

# jawinski_markett
#brainage2025.full.eur.gm.gz
#brainage2025.full.eur.gwm.gz
#brainage2025.full.eur.wm.gz
python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/jawinski_markett/brainage2025.full.eur.gm.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/jawinski_markett_gm_gwas_harmonised.txt.gz \
	-separator " " \
	--chromosome_format \
	-output_column_map CHR chromosome \
	-output_column_map BP position \
	-output_column_map ID variant_id \
	-output_column_map A1 effect_allele \
	-output_column_map A2 non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 54890 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/jawinski_markett/brainage2025.full.eur.wm.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/jawinski_markett_wm_gwas_harmonised.txt.gz \
	-separator " " \
	--chromosome_format \
	-output_column_map CHR chromosome \
	-output_column_map BP position \
	-output_column_map ID variant_id \
	-output_column_map A1 effect_allele \
	-output_column_map A2 non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 54890 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pvalue zscore effect_size standard_error sample_size n_cases frequency

python $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/gwas_data/jawinski_markett/brainage2025.full.eur.gwm.gz \
	-liftover $LIFTOVER/hg19ToHg38.over.chain.gz \
	-snp_reference_metadata $DATA/gtex_reference/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz METADATA \
	-output $OUTPUT/jawinski_markett_gwm_gwas_harmonised.txt.gz \
	-separator " " \
	--chromosome_format \
	-output_column_map CHR chromosome \
	-output_column_map BP position \
	-output_column_map ID variant_id \
	-output_column_map A1 effect_allele \
	-output_column_map A2 non_effect_allele \
	-output_column_map BETA effect_size \
	-output_column_map SE standard_error \
	-output_column_map P pvalue \
	--insert_value sample_size 54890 \
	-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele pvalue zscore effect_size standard_error sample_size n_cases frequency


##########################################################################
# IMPUTATION
for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/pineiro_fjell_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/pineiro_fjell_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/kim_lee_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/kim_lee_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/zhao_zhao_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/zhao_zhao_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/kaufmann_westlye_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/kaufmann_westlye_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/fan_huang_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/fan_huang_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/leonardsen_wang_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/leonardsen_wang_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/wen_davatzikos_wm_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/wen_davatzikos_wm_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/wen_davatzikos_gm_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/wen_davatzikos_gm_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/wen_davatzikos_fc_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/wen_davatzikos_fc_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/wen_davatzikos_combined_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/wen_davatzikos_combined_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

for chr in {1..22}
do
	for sub_batch in {0..9}
	do
		echo "~~~~~~~~~~ imputation chromosome ${chr} sub-batch ${sub_batch} ~~~~~~~~~~"
		python $GWAS_TOOLS/gwas_summary_imputation.py \
			-by_region_file $DATA/eur_ld.bed.gz \
			-gwas_file $OUTPUT/jawinski_markett_gwm_gwas_harmonised.txt.gz \
			-parquet_genotype $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.chr${chr}.variants.parquet \
			-parquet_genotype_metadata $ROOT/genotype/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.variants_metadata.parquet \
			-window 100000 \
			-parsimony 7 \
			-chromosome ${chr} \
			-regularization 0.1 \
			-frequency_filter 0.01 \
			-sub_batches 10 \
			-sub_batch ${sub_batch} \
			--standardise_dosages \
			--cache_variants \
			-output $OUTPUT/jawinski_markett_gwm_gwas_harmonised_imputed_chr${chr}_sb${sub_batch}_reg0.1_ff0.01_by_region.txt.gz
	done
done

#########################################################################################
# merge imputed sumstats batches into one file

python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
	-gwas_file $OUTPUT/pineiro_fjell_gwas_harmonised.txt.gz \
	-folder $OUTPUT \
	-pattern pineiro_fjell_gwas_harmonised_imputed_.* \
	-parsimony 7 \
	-output $OUTPUT/pineiro_fjell_gwas_harmonised_imputed_merged.txt.gz

############################################################################################
#eqtl
#	individual
#	brain multi
#	all multi
#sqtl
#	individual
#	brain multi
#	all multi

# Run S-PrediXcan for all tissues individually
study="pineiro_fjell"
n_size=35712
MODEL="$DATA/models/sqtl/mashr"
for this_model in `ls $MODEL/*.db -1 | xargs -n 1 basename | sed 's/\.[^.]*$//'`
do
	echo "model: $this_model"
	python $METAXCAN/SPrediXcan.py \
		--gwas_file $OUTPUT/${study}_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path $MODEL/${this_model}.db \
		--covariance $MODEL/${this_model}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--throw \
		--gwas_N $n_size \
		--chromosome_column chromosome \
		--position_column position \
		--freq_column frequency \
		--zscore_column zscore \
		--pvalue_column pvalue \
		--output_file $OUTPUT/sqtl/sqtl_spredixcan_${study}_${this_model}.csv
done

# run s metaxcan for all
# expression - all tissues
MODEL="$DATA/models/eqtl/mashr"
for study in yi_huang jawinski_markett_gm jawinski_markett_wm jawinski_markett_gwm kaufmann_westlye zhao_zhao leonardsen_wang kim_lee pineiro_fjell
do
	echo $study
	python $METAXCAN/SMulTiXcan.py \
		--models_folder $MODEL \
		--models_name_pattern "mashr_(.*).db" \
		--snp_covariance $DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
		--metaxcan_folder $OUTPUT/ \
		--metaxcan_filter "spredixcan_${study}_(.*).csv" \
		--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
		--gwas_file $OUTPUT/${study}_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--keep_non_rsid \
		--model_db_snp_key varID \
		--cutoff_condition_number 30 \
		--verbosity 7 \
		--throw \
		--output $OUTPUT/smultixcan_${study}_all_tissues.csv
done

# splicing - all tissues
MODEL="$DATA/models/sqtl/mashr"
for study in yi_huang jawinski_markett_gm jawinski_markett_wm jawinski_markett_gwm kaufmann_westlye zhao_zhao leonardsen_wang kim_lee pineiro_fjell
do
	echo $study
	python $METAXCAN/SMulTiXcan.py \
		--models_folder $MODEL \
		--models_name_pattern "mashr_(.*).db" \
		--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
		--metaxcan_folder $OUTPUT/sqtl/ \
		--metaxcan_filter "sqtl_spredixcan_${study}_(.*).csv" \
		--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
		--gwas_file $OUTPUT/${study}_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--keep_non_rsid \
		--model_db_snp_key varID \
		--cutoff_condition_number 30 \
		--verbosity 7 \
		--throw \
		--output $OUTPUT/sqtl/sqtl_smultixcan_${study}_all_tissues.csv
done

# expression - brain tissues
MODEL="$DATA/models/eqtl/mashr"
for study in yi_huang jawinski_markett_gm jawinski_markett_wm jawinski_markett_gwm kaufmann_westlye zhao_zhao leonardsen_wang kim_lee pineiro_fjell
do
	echo $study
	python $METAXCAN/SMulTiXcan.py \
		--models_folder $MODEL \
		--models_name_filter "mashr_Brain_(.*).db" \
		--snp_covariance $DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
		--metaxcan_folder $OUTPUT/ \
		--metaxcan_filter "spredixcan_${study}_mashr_Brain_(.*).csv" \
		--metaxcan_file_name_parse_pattern "(.*)_${study}_(.*).csv" \
		--gwas_file $OUTPUT/${study}_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--keep_non_rsid \
		--model_db_snp_key varID \
		--cutoff_condition_number 30 \
		--verbosity 7 \
		--throw \
		--output $OUTPUT/smultixcan_${study}_brain_tissues.csv
done

# splicing - brain tissues
MODEL="$DATA/models/sqtl/mashr"
for study in yi_huang jawinski_markett_gm jawinski_markett_wm jawinski_markett_gwm kaufmann_westlye zhao_zhao leonardsen_wang kim_lee pineiro_fjell
do
	echo $study
	python $METAXCAN/SMulTiXcan.py \
		--models_folder $MODEL \
		--models_name_filter "mashr_Brain_(.*).db" \
		--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
		--metaxcan_folder $OUTPUT/sqtl/ \
		--metaxcan_filter "sqtl_spredixcan_${study}_mashr_Brain_(.*).csv" \
		--metaxcan_file_name_parse_pattern "(.*)_${study}_(.*).csv" \
		--gwas_file $OUTPUT/${study}_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--keep_non_rsid \
		--model_db_snp_key varID \
		--cutoff_condition_number 30 \
		--verbosity 7 \
		--throw \
		--output $OUTPUT/sqtl/sqtl_smultixcan_${study}_brain_tissues.csv
done

# smultixcan brain
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_filter "mashr_Brain_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/ \
	--metaxcan_filter "spredixcan_wen_davatzikos_fc_mashr_Brain_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_fc_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_fc_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/smultixcan_wen_davatzikos_fc_brain_tissues.csv

############### sQTLs
MODEL="$DATA/models/sqtl/mashr"
# all tissues individually
for this_model in `ls $MODEL/*.db -1 | xargs -n 1 basename | sed 's/\.[^.]*$//'`
do
	echo "model: $this_model"
	python $METAXCAN/SPrediXcan.py \
		--gwas_file $OUTPUT/wen_davatzikos_gm_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path $MODEL/${this_model}.db \
		--covariance $MODEL/${this_model}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--throw \
		--gwas_N 31129 \
		--chromosome_column chromosome \
		--position_column position \
		--freq_column frequency \
		--zscore_column zscore \
		--pvalue_column pvalue \
		--output_file $OUTPUT/sqtl/sqtl_spredixcan_wen_davatzikos_gm_${this_model}.csv
		
		python $METAXCAN/SPrediXcan.py \
		--gwas_file $OUTPUT/wen_davatzikos_wm_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path $MODEL/${this_model}.db \
		--covariance $MODEL/${this_model}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--throw \
		--gwas_N 31129 \
		--chromosome_column chromosome \
		--position_column position \
		--freq_column frequency \
		--zscore_column zscore \
		--pvalue_column pvalue \
		--output_file $OUTPUT/sqtl/sqtl_spredixcan_wen_davatzikos_wm_${this_model}.csv
		
		python $METAXCAN/SPrediXcan.py \
		--gwas_file $OUTPUT/wen_davatzikos_fc_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path $MODEL/${this_model}.db \
		--covariance $MODEL/${this_model}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--throw \
		--gwas_N 31129 \
		--chromosome_column chromosome \
		--position_column position \
		--freq_column frequency \
		--zscore_column zscore \
		--pvalue_column pvalue \
		--output_file $OUTPUT/sqtl/sqtl_spredixcan_wen_davatzikos_fc_${this_model}.csv
		
		python $METAXCAN/SPrediXcan.py \
		--gwas_file $OUTPUT/wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz \
		--snp_column panel_variant_id \
		--effect_allele_column effect_allele \
		--non_effect_allele_column non_effect_allele \
		--zscore_column zscore \
		--model_db_path $MODEL/${this_model}.db \
		--covariance $MODEL/${this_model}.txt.gz \
		--keep_non_rsid \
		--additional_output \
		--model_db_snp_key varID \
		--throw \
		--gwas_N 31129 \
		--chromosome_column chromosome \
		--position_column position \
		--freq_column frequency \
		--zscore_column zscore \
		--pvalue_column pvalue \
		--output_file $OUTPUT/sqtl/sqtl_spredixcan_wen_davatzikos_combined_${this_model}.csv
done

############################################################################################
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_pattern "mashr_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_gm_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_gm_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_gm_all_tissues.csv

python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_pattern "mashr_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_wm_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_wm_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_wm_all_tissues.csv
	
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_pattern "mashr_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_fc_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_fc_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_fc_all_tissues.csv
	
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_pattern "mashr_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_combined_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_combined_all_tissues.csv

# smultixcan brain 
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_filter "mashr_Brain_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_wm_mashr_Brain_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_wm_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_wm_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_wm_brain_tissues.csv
	
python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_filter "mashr_Brain_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_gm_mashr_Brain_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_gm_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_gm_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_gm_brain_tissues.csv	

python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_filter "mashr_Brain_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_fc_mashr_Brain_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_fc_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_fc_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_fc_brain_tissues.csv

python $METAXCAN/SMulTiXcan.py \
	--models_folder $MODEL \
	--models_name_filter "mashr_Brain_(.*).db" \
	--snp_covariance $DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz \
	--metaxcan_folder $OUTPUT/sqtl/ \
	--metaxcan_filter "sqtl_spredixcan_wen_davatzikos_combined_mashr_Brain_(.*).csv" \
	--metaxcan_file_name_parse_pattern "(.*)_combined_(.*).csv" \
	--gwas_file $OUTPUT/wen_davatzikos_combined_gwas_harmonised_imputed_merged.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--zscore_column zscore \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--output $OUTPUT/sqtl/sqtl_smultixcan_wen_davatzikos_combined_brain_tissues.csv

