source params.cfg

# TODO: At this moment you should check the scaffolds manually and save the list of accepted scaffolds in spades_scaffolds_selected.txt
xargs seqkit faidx -f 05_spades_final/scaffolds.fasta < spades_scaffolds_selected.txt > 05_spades_final/scaffolds_selected.fasta
quast.py -1 04_filtered/PoV-01B_R1_cutadapt_filtered.fq.gz -2 04_filtered/PoV-01B_R2_cutadapt_filtered.fq.gz -o 05_spades_final/quast 05_spades_final/scaffolds_selected.fasta

prokka --outdir output --evalue 1e-5 --coverage 60 --force --gcode 1 --addgenes --prefix PoV-01B.prokka --compliant --kingdom Viruses --genus 'Pyramimonas orientalis' --species virus --strain 01B --locustag POV01B 05_spades_final/scaffolds_selected.fasta

perl scripts/biotags.pl -i 06_annotate/PoV-01B.prokka.gbk -T id,seq | awk '{$1 = sprintf("%s Pyramimonas orientalis virus 01B genomic scaffold %s", $1, $1)}1' OFS=\\t | \
	seqkit tab2fx > 07_submit/PoV-01B.submission.fsa

gmsn.pl --format GFF3 --virus --output 06_annotate/PoV-01B.gmsn.gff3 07_submit/PoV-01B.submission.fsa
yes O | tRNAscan-SE -o 06_annotate/PoV-01B.trna -f 06_annotate/PoV-01B.trna.struct -G 07_submit/PoV-01B.submission.fsa

perl scripts/annotate-prokka.pl -p 06_annotate/PoV-01B.prokka.gbk -g 06_annotate/PoV-01B.gmsn.gff3 -t 06_annotate/PoV-01B.trna -m 06_annotate/PoV-01B.manual.gbk -l "$locus_tag" > 06_annotate/PoV-01B.merged.gbk

perl scripts/gbk2tbl.pl < 06_annotate/PoV-01B.merged.gbk > 07_submit/PoV-01B.submission.tbl

tbl2asn -t 07_submit/PoV-01B.submission.sbt -a s -i 07_submit/PoV-01B.submission.fsa -V vb

ln -fs ../07_submit/PoV-01B.submission.fsa output/PoV-01B.fna
ln -fs ../07_submit/PoV-01B.submission.gbf output/PoV-01B.gbk

perl scripts/biotags.pl -i output/PoV-01B.gbk -p CDS -t locus_tag,translation | seqkit tab2fx > output/PoV-01B.faa

awk '$1=="LOCUS"{f=sprintf("output/split/%s.gbk",$2)}{print > f}' output/PoV-01B.gbk
