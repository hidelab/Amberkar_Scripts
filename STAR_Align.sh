tdp_r1_path="/fastdata/md4zsa/TDP_Omics_Study/Data/R1"
tdp_r2_path="/fastdata/md4zsa/TDP_Omics_Study/Data/R2"
tdp_r1_files=($(find $tdp_r1_path -type f|sort -n|awk '{print $0,"\n"}'))
tdp_r2_files=($(find $tdp_r2_path -type f|sort -n|awk '{print $0,"\n"}'))
out_dir="/fastdata/md4zsa/TDP_Omics_Study/Results/STAR_Aligner"
star_index="/fastdata/md4zsa/Genome_Indices/STAR_Index/"
hs_ann="/fastdata/md4zsa/Ensembl_Genomes/annotation/Homo_sapiens.GRCh38.82.gtf"
genome_dir="/fastdata/md4zsa/Genome_Indices/STAR_Index/"
for i in `seq 0 39`;
do
STAR --runThreadN 32 --genomeDir $genome_dir --readFilesIn ${tdp_r1_files[i]} ${tdp_r2_files[i]}--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $out_dir/$(basename ${tdp_r1_files[i]}|cut -d'_' -f1) --sjdbGTFfile $hs_ann
done