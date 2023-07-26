


#############################################
# RNASEQ ANALYSIS FOR CANDIDA PARAPSILOSIS  #
#############################################



# Introduction: Data for this analysis was downloaded from SRA.
# The experiment aims to study gene expression in Candida parapsilosis under two growth conditions: planktonic and biofilm. There are 6 fastq files, corresponding to 3 replicates of each treatment.
# Accession numbers for the files:"SRR1278968" "SRR1278969" "SRR1278970" "SRR1278971" "SRR1278973" "SRR1278972".


#############################################
# DOWNLOAD FILES                            #
#############################################


# To download and split the files one at a time:
# Example:
# fastq-dump SRR1278968 --split-files

# To download all the files simultaneously, create a .sh file with the following code:

#!/bin/bash

# Array of replicas
replicas=("SRR1278968" "SRR1278969" "SRR1278970" "SRR1278971" "SRR1278973" "SRR1278972")

for replica in "${replicas[@]}"; do
  # Create the full path to the replica folder
  carpeta_replica="${replica}/trim"
  
  # Create the folder if it does not exist
  mkdir -p "${carpeta_replica}"
  
  # Execute the fastq-dump command
  fastq-dump "${replica}" --split-files
done

# Change the permissions of the file
# chmod +x archivo.sh

# Execute the file
# ./archivo.sh



#############################################
# QUALITY CONTROL CHECK                     #
#############################################

# The quality control can also be done for each file one at a time.
# Example:
# fastqc SRR1278968_1.fastq
# fastqc SRR1278968_2.fastq

# To perform quality control for all files simultaneously create a file:

#!/bin/bash

# Array of directories for the folders

carpetas=(
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278969"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278970"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278971"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278972"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278973"
)

# Iterate over the folders
for carpeta in "${carpetas[@]}"; do
  # Get the folder name
  nombre_carpeta=$(basename "$carpeta")
  
  # Paths to the files
  archivo_1="$carpeta/${nombre_carpeta}_1.fastq"
  archivo_2="$carpeta/${nombre_carpeta}_2.fastq"
  
  # Execute fastqc
  fastqc "$archivo_1"
  fastqc "$archivo_2"
done




#############################################
# CLEANING                                  #
#############################################

# The cleaning can also be done for each sequence one at a time.
# Example:
# skewer -m pe -l 36 -q 15 -Q 15 -o trim/SRR1278968 -t 4 SRR1278968_1.fastq SRR1278968_2.fastq
# fastqc SRR1278968-trimmed-pair1.fastq
# fastqc SRR1278968-trimmed-pair2.fastq

# To perform cleaning for all sequences simultaneously:

#!/bin/bash

# Array of directories for the folders
carpetas=(
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278969"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278970"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278971"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278972"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278973"
)

# Iterate over the folders
for carpeta in "${carpetas[@]}"; do
  # Get the folder name
  nombre_carpeta=$(basename "$carpeta")
  
  # Paths to the files
  archivo_1="$carpeta/${nombre_carpeta}-trimmed-pair1.fastq"
  archivo_2="$carpeta/${nombre_carpeta}-trimmed-pair2.fastq"
  
  # Create an output directory for skewer
  carpeta_trim="trim/$nombre_carpeta"
  mkdir -p "$carpeta_trim"

  # Execute skewer
  skewer -m pe -l 36 -q 15 -Q 15 -o "$carpeta_trim" -t 4 "$archivo_1" "$archivo_2"
done

# Quality check for the cleaned fastq files
for carpeta in "${carpetas[@]}"; do
  # Get the folder name
  nombre_carpeta=$(basename "$carpeta")
  
  # Paths to the files
  archivo_1="$carpeta/${nombre_carpeta}-trimmed-pair1.fastq"
  archivo_2="$carpeta/${nombre_carpeta}-trimmed-pair2.fastq"
  
  # Execute fastqc
  fastqc "$archivo_1"
  fastqc "$archivo_2"
done


#############################################
# DOWNLOADING AND INDEXING OF REFERENCE GENOME AND GTF FILE #
#############################################

# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/5480/

# Genome index generation using STAR software
#(“splicing-aware” aligner that can recognize the difference between a read aligning across an exon–intron boundary and a read with a short insertion.)

STAR --runMode genomeGenerate --genomeDir STAR_genome_gencode --genomeFastaFiles GCA_000182765.2_ASM18276v2_genomic.fna --sjdbGTFfile genomic.gtf --genomeSAindexNbases 10



#############################################
# ALIGNMENT WITH STAR SOFTWARE              #
#############################################

# For a single file
# Example:
# STAR --genomeDir STAR_genome_gencode --readFilesIn SRR1278968-trimmed-pair1.fastq SRR1278968-trimmed-pair2.fastq --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 10 --limitBAMsortRAM 1152904740

# For all files

#!/bin/bash

# Array of directories for the folders
carpetas=(
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278969"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278970"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278971"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278972"
  "/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/SRR1278973"
)

# Iterate over the folders
for carpeta in "${carpetas[@]}"; do
  # Get the folder name
  nombre_carpeta=$(basename "$carpeta")

  # Directory for STAR output
  output_dir="$carpeta/STAR_output"
  mkdir -p "$output_dir"
  
  # Execute STAR
  genomeDir="STAR_genome_gencode"
  archivo_1="$carpeta/${nombre_carpeta}_1.fastq"
  archivo_2="$carpeta/${nombre_carpeta}_2.fastq"
  
  STAR --genomeDir "$genomeDir" --readFilesIn "$archivo_1" "$archivo_2" --outSAMtype BAM SortedByCoordinate --genomeSAindexNbases 10 --limitBAMsortRAM 1152904740 --outFileNamePrefix "$output_dir/${nombre_carpeta}_"
done


#############################################
# TABLES OF GENE EXPRESSION COUNTS          #
#############################################

# For a single file
# Example:
# htseq-count -f bam -r pos -t exon -i gene_id -m intersection-strict Aligned.sortedByCoord.out.bam GCF_000182765.1_ASM18276v2_genomic.gtf > counts.txt

# Concatenate count tables

# Create a script named merge_counts.sh
nano merge_counts.sh

#!/bin/bash

# List of input files to merge
input_files=("SRR1278968_counts.txt" "SRR1278969_counts.txt" "SRR1278970_counts.txt" "SRR1278971_counts.txt" "SRR1278972_counts.txt" "SRR1278973_counts.txt")

# Output file name for merged data
output_file="merged_counts.txt"

# Sort and merge the files
sort -k1,1 "${input_files[0]}" > "$output_file"

for ((i = 1; i < ${#input_files[@]}; i++)); do
  join -a 1 -a 2 -e "0" -o auto "$output_file" <(sort -k1,1 "${input_files[i]}") > merged_tmp.txt
  mv merged_tmp.txt "$output_file"
done

chmod +x merge_counts.sh

./merge_counts.sh






