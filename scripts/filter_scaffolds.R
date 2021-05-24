library(dplyr)
library(tidyr)
library(bioformatr)
library(Biostrings)
library(IRanges)
library(mmgenome2)

options(stringsAsFactors = F)

args <- commandArgs(T)

scaffolds.fname <- args[1]
viruses.fname   <- args[2]
sanger.fname    <- args[3]
bacteria.fname  <- args[4]

# scaffolds.fname <- "04_filtered/scaffolds.fasta"
# viruses.fname   <- "04_filtered/scaffolds_mesomimivirinae.blast"
# sanger.fname    <- "04_filtered/scaffolds_sanger.blast"
# bacteria.fname  <- "04_filtered/scaffolds_contamination.blast"

viruses <- read.outfmt6(viruses.fname, col.names.pre = "OG") %>%
	distinct(OG, sseqid) %>%
	group_by(sseqid) %>%
	summarize(OGs = n())

sanger <- read.outfmt6(sanger.fname) %>%
	filter(pident > 97) %>%
	select(qseqid, sseqid, pident, length)

bacteria <- read.outfmt6(bacteria.fname) %>%
	split(.$qseqid) %>%
	lapply(function(x) {
		qstart <- ifelse(x$qstart < x$qend, x$qstart, x$qend)
		qend   <- ifelse(x$qstart > x$qend, x$qstart, x$qend)
		IRanges(start = qstart, end = qend) %>%
			coverage %>%
			{crossprod(.@lengths, .@values > 0)}
	}) %>%
	{data.frame(qseqid = names(.), bact.cov = unlist(.))}

all.scaffolds <- readDNAStringSet(scaffolds.fname) %>%
	`[`(width(.) >= 500)

scaffolds <- data.frame(qseqid = names(all.scaffolds)) %>%
	extract(qseqid, into = c("num", "scaffold.len", "cov"), regex = "NODE_([0-9]+)_length_([0-9]+)_cov_([0-9.]+)", convert = T, remove = F) %>%
	left_join(bacteria, by = "qseqid") %>%
	replace_na(list(bact.cov = 0)) %>%
	mutate(bact.pct = bact.cov / scaffold.len) %>%
	left_join(sanger, by = "qseqid") %>%
	left_join(viruses, by = c(qseqid = "sseqid"))

sc.taxonomy <- scaffolds %>%
	mutate(is.bact = bact.pct > 0.1 | bact.cov > 1000, is.vir = bact.cov == 0 & OGs > 3, is.sanger = !is.na(sseqid)) %>%
	filter(is.bact | is.vir) %>%
	mutate(tax = ifelse(is.bact, "Bacteria", "Viruses")) %>%
	distinct(qseqid, tax)

sc.sanger <- distinct(sanger, qseqid) %>% pull

sc.coverage <- distinct(scaffolds, qseqid, cov)

sc.sanger <- mutate(scaffolds, sanger = ifelse(is.na(sseqid), NA, "Sanger")) %>% distinct(qseqid, sanger)

mm <- mmload(
	assembly = all.scaffolds,
	coverage = sc.coverage,
	taxonomy = sc.taxonomy,
	additional = sc.sanger,
	verbose = T,
	kmer_pca = T,
	kmer_BH_tSNE = F
)

if (interactive()) {
	# if working interactively - create the selection area and save as a data frame in filter_scaffolds_1.txt
	mmplot(mm.sel, x = "cov_cov", y = "gc", color_by = "tax", locator = T)
}

selection <- read.table("filter_scaffolds_1.txt", col.names = c("cov_cov", "gc"))

g <- mmplot(mm, x = "cov_cov", y = "gc", color_by = "tax", x_scale = "log10", locator = F, selection = selection)
ggsave("cov_gc.svg", g)

mm_subset <- mmextract(mm, selection = selection) %>% filter(is.na(tax) | tax != "Bacteria")

selected.scaffolds <- all.scaffolds[mm_subset$scaffold]

rm(assembly)
mm.sel <- mmload(
	assembly = selected.scaffolds,
	coverage = sc.coverage,
	taxonomy = sc.taxonomy,
	additional = sc.sanger,
	verbose = T,
	kmer_pca = T,
	kmer_BH_tSNE = F
)

if (interactive()) {
	# if working interactively - create the selection area and save as a data frame in filter_scaffolds_2.txt
	mmplot(mm.sel, x = "PC1", y = "PC2", color_by = "sanger", locator = T)
}

selection <- read.table("filter_scaffolds_2.txt", col.names = c("PC1", "PC2"))

g <- mmplot(mm.sel, x = "PC1", y = "PC2", color_by = "sanger", locator = F, selection = selection)
ggsave("PC1_PC2.svg", g)

mmextract(mm.sel, selection = selection) %>% pull(scaffold) %>% cat(sep = "\n")
