# Pyramimonas orientalis virus genome, assembly pipeline

This is the pipeline for the assembly of the Pyramimonas orientalis virus genome. The raw data are available from [PRJNA641252](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA641252). See the paper for [details](http://dx.doi.org/10.1016/j.cub.2020.09.056).

## Dependencies

* R, perl
* [cutadapt](https://cutadapt.readthedocs.io/)
* [pilon](https://github.com/broadinstitute/pilon) - the path `$PILON_path` should be specified in params.cfg
* [prokka](https://github.com/tseemann/prokka)
* [NCBI blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [spades](https://github.com/ablab/spades)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [samtools](http://www.htslib.org/)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
* [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)
* [bioformatr](https://github.com/alephreish/bioformatr/)
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
* [mmgenome2](https://github.com/KasperSkytte/mmgenome2)
* [Bio::SeqIO](https://metacpan.org/pod/Bio::SeqIO)
* [Bio::Tools::GFF](https://metacpan.org/pod/Bio::Tools::GFF)

The blast+ database of [SILVA\_\*\_SSURef\_tax\_silva.fasta](https://www.arb-silva.de/no_cache/download/archive/) should be availbale from the path `$SILVA_path` specified in params.cfg.

## Usage

The script `assemble.sh` assembles the data. Whenever data changes, the files `filter_scaffolds_1.txt` and `filter_scaffolds_2.txt` should be changed using the interactive mode of `mmgenome2` (refer to `filter_scaffolds.R` for details).

The script `finalize.sh` produces the final output files and relies on the user-supplied list of selected scaffolds `spades_scaffolds_selected.txt`.
