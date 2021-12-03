#!/bin/bash

#################################################
################ Load variables  ################
#################################################

config=$1
. $config

dir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo $dir
echo $n
echo $typeS
echo $output

mkdir $output
###############################################
############# Load env and module #############
###############################################
echo 'activate TRAWLING_env...'
conda activate TRAWLING_env 

########################################################
####### Download gtf and fasta file from gencode #######
########################################################
echo 'download gtf and fasta file from gencode...'

mkdir $output/'reference/'
GTF_file=$output/'reference/gencode.gtf' 
wget -P $output/'reference' https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz -O $GTF_file'.gz'
gunzip $GTF_file'.gz'
Ref_file=$output/'reference/homo-sapiens-gencode.fa' 
wget -P $output/'reference' https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz -O $Ref_file'.gz'
gunzip $Ref_file'.gz'

#####################################################
############# Create gtfs and beds file #############
#####################################################
subdir_tmp=$output/'tmp/'
mkdir $subdir_tmp
echo 'Create gtfs and beds file...'
R CMD BATCH --vanilla "--args $GTF_file $subdir_tmp $n $typeS" $dir/script/create_gtf_bed_files.R

################################################################
############# Extract donor/acceptor splice motifs #############
################################################################
echo 'Extract donor/acceptor splice motifs...'
mkdir $output'/tmp/fasta_bed/'
bedtools getfasta -name -s -fi $Ref_file -bed $subdir_tmp'/bed_tmp/homo_sapiens_gencode_exons.bed' > $subdir_tmp'/fasta_bed/homo_sapiens_gencode_exons.fa'
awk 'BEGIN{RS=">";OFS="";}NR>1{print $1 "\t" $2}' $subdir_tmp'/fasta_bed/homo_sapiens_gencode_exons.fa' |  awk -F' ' ' {{split($1, arr, "::"); print arr[1] "\t" $2}}' |  awk '{n=split($1,arr,"."); pfx=""; for (i=1;i<n;i++) {printf "%s%s",pfx,arr[i]; pfx="."}; print "\t" arr[n] "\t" $2}' >  $subdir_tmp'/fasta_bed/homo_sapiens_gencode_exons2.fa'


bedtools getfasta -name -s -fi  $Ref_file -bed $subdir_tmp'/bed_tmp/homo_sapiens_gencode_jnc.bed' > $subdir_tmp'/fasta_bed/homo_sapiens_gencode_jnc.fa'
awk 'BEGIN{RS=">";OFS="";}NR>1{print $1 "\t" $2}' $subdir_tmp'/fasta_bed/homo_sapiens_gencode_jnc.fa' |  awk -F' ' ' {{split($1, arr, "::"); print arr[1] "\t" $2}}' |  awk '{n=split($1,arr,"."); pfx=""; for (i=1;i<n;i++) {printf "%s%s",pfx,arr[i]; pfx="."}; print "\t" arr[n] "\t" $2}' > $subdir_tmp'/fasta_bed/homo_sapiens_gencode_jnc2.fa'


######################################################
############# Merge the information and  #############
############# create TRAWLING's gtf file #############
######################################################
echo 'Create TRAWLING gtf file...'
R CMD BATCH --vanilla "--args $subdir_tmp $output $typeS" $dir/script/create_TRAWLING_gtf.R

######################################################
################# Extract fasta file #################
######################################################
echo 'Extract sequences fasta format...'
TRAWLING_gtf=$output'/gtf/TRAWLING.gtf' 
TRAWLING_fasta=$output'/fasta/TRAWLING-cdna'$typeS'.fa'
Gencode_gtf=$output'/gtf/HS_gencode.gtf' 
Gencode_fasta=$output'/fasta/HS_gencode-cdna'$typeS'.fa'

gffread -w $TRAWLING_fasta -g $Ref_file $TRAWLING_gtf
gffread -w $Gencode_fasta -g $Ref_file $Gencode_gtf

