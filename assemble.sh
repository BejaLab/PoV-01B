
source params.cfg

mkdir -p 02_data 04_filtered output

# Prepare data
cat 01_miseq/PoV-01B_S*_R1_001.fastq.gz > 02_data/PoV-01B_R1.fq.gz
cat 01_miseq/PoV-01B_S*_R2_001.fastq.gz > 02_data/PoV-01B_R2.fq.gz

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o 02_data/PoV-01B_R1_cutadapt.fq.gz -p 02_data/PoV-01B_R2_cutadapt.fq.gz 02_data/PoV-01B_{R1,R2}.fq.gz &> 02_data/cutadapt.log

# Prepare Sanger
seqkit fx2tab 01_sanger/phrap/edit_dir/PoV.contigs | cut -f2 | tr x n | nl -w1 | sed s/^/CONTIG_/ | seqkit tab2fx > 01_sanger/phrap.fasta
makeblastdb -in 01_sanger/phrap.fasta -dbtype nucl
bowtie2-build 01_sanger/phrap.fasta{,}
bowtie2 --very-sensitive --no-unal -x 01_sanger/phrap.fasta -1 02_data/PoV-01B_R1_cutadapt.fq.gz -2 02_data/PoV-01B_R2_cutadapt.fq.gz | samtools view -Sb | samtools sort > 01_sanger/phrap.bam
samtools index 01_sanger/phrap.bam
java -jar "$PILON_path/pilon.jar" --frags 01_sanger/phrap.bam --genome 01_sanger/phrap.fasta --outdir 02_data --output phrap_pilon

# Do meta spades
spades.py -1 02_data/PoV-01B_R1_cutadapt.fq.gz -2 02_data/PoV-01B_R2_cutadapt.fq.gz --meta -o 03_spades_meta -t "$CPUs"

# Do blasts against the MG assembly
makeblastdb -in 03_spades_meta/scaffolds.fasta -dbtype nucl
blastn -num_threads 20 -max_target_seqs 1 -query 03_spades_meta/scaffolds.fasta -db "$SILVA_path" -evalue 1e-10 -outfmt 6 > 04_filtered/scaffolds_SILVA.blast
parallel --tagstring {/.} tblastn -num_threads "$CPUs" -db 03_spades_meta/scaffolds.fasta -query {} -outfmt 6 -evalue 1e-5 ::: 00_refs/mesomimivirinae/*.faa > 04_filtered/scaffolds_mesomimivirinae.blast
blastn -max_target_seqs 1 -num_threads "$CPUs" -query 03_spades_meta/scaffolds.fasta -db 00_refs/bacteria/contamination.fasta -outfmt 6 -evalue 1e-10 > 04_filtered/scaffolds_contamination.blast
blastn -num_threads "$CPUs" -db 01_sanger/phrap.fasta -query 03_spades_meta/scaffolds.fasta -outfmt 6 -evalue 1e-20 > 04_filtered/scaffolds_sanger.blast

# Filter the scaffolds
# TODO: filter_scaffolds.R might need user intervention in interactive mode
Rscript scripts/filter_scaffolds.R 03_spades_meta/scaffolds.fasta 04_filtered/scaffolds_{mesomimivirinae,sanger,contamination}.blast | xargs samtools faidx 03_spades_meta/scaffolds.fasta > 04_filtered/scaffolds_filtered.fasta

# Recruit reads against the selected scaffolds
bowtie2-build 04_filtered/scaffolds_filtered.fasta{,}
bowtie2 -p "$CPUs" --very-sensitive-local --no-unal -x 04_filtered/scaffolds_filtered.fasta -1 02_data/PoV-01B_R1_cutadapt.fq.gz -2 02_data/PoV-01B_R2_cutadapt.fq.gz -S 04_filtered/scaffolds_filtered.sam
cut -f1 04_filtered/scaffolds_filtered.sam | grep -v ^@ | uniq | seqkit grep -f- 02_data/PoV-01B_R1_cutadapt.fq.gz | gzip > 04_filtered/PoV-01B_R1_cutadapt_filtered.fq.gz
cut -f1 04_filtered/scaffolds_filtered.sam | grep -v ^@ | uniq | seqkit grep -f- 02_data/PoV-01B_R2_cutadapt.fq.gz | gzip > 04_filtered/PoV-01B_R2_cutadapt_filtered.fq.gz

# Do hybrid spades
spades.py -k 21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97,101,105,109,113,117,121,125 --trusted-contigs 02_data/phrap_pilon.fasta -1 04_filtered/PoV-01B_R1_cutadapt_filtered.fq.gz -2 04_filtered/PoV-01B_R2_cutadapt_filtered.fq.gz -o spades_final -t "$CPUs"
