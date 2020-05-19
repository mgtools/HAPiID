#!/bin/bash
#
#Script where I will get all the potential ORFs with peptide evidence, deduplicate them using CD-hit at 100%
#and then will extract the corresponding contigs that they come from, predict the proteins from those 
#contigs using various predictors (FGS, prodigal, genemark) and construct a blast protein database
#using the three predicted source and search these potential ORFs against them
#

START=$(date +%s)


genomes_dir='../data/genomes/'
contigs2genomes_f='../data/contigs2genomes.json'
proteomes_dir='../data/proteomes/'

while getopts i:o:t:n:h: option
do
case "${option}"
in
    i) expressedORFs_dir=${OPTARG};;
    o) out_dir=${OPTARG};;
    t) n_cores=${OPTARG};;
    n) n_topGenomes=${OPTARG};;
    h) HAPiID_n_topGenomes=${OPTARG};;
esac
done


fgs_out=$out_dir"FGS_out/"
prodigal_out=$out_dir"prodigal_out/"
genemark_out=$out_dir"genemark_out/"

blast_out=$out_dir"blast_out/"



if [ ! -d $out_dir ]; then
    mkdir $out_dir
fi

#concatinate all the potential ORFs with protein evidence into one file
#cat $expressedORFs_dir*'_top'$HAPiID_n_topGenomes'_top'$n_topGenomes'_'*.fasta > $out_dir"allExpressedORFsSoFar.fasta"
#deduplicate the potential ORFs by cd-hit at 100%
#cd-hit -i $out_dir"allExpressedORFsSoFar.fasta" -o $out_dir"allExpressedORFsSoFar_c1.0.fasta" -n 5 -M 16000 -d 0 -T $n_cores
#extract the contigs which these ORFs are present on
#python3 get_ORFcontigs.py $out_dir"allExpressedORFsSoFar_c1.0.fasta" $out_dir
#now predicit protein coding genes using different predictors
if [ ! -d $fgs_out ]; then
    mkdir $fgs_out 
fi

if [ ! -d $prodigal_out ]; then
    mkdir $prodigal_out
fi

if [ ! -d $genemark_out ]; then
    mkdir $genemark_out
fi

#genmark predict
echo "running genemark.hmm"
#gmhmmp -m ~/tree_match_pipeline/tree_match_prereqs/genemark_suite_linux_64/gmsuite/heuristic_mod/heu_11_47.mod $out_dir"contigsWithPeptides.fasta" -o $genemark_out"contigsWithPeptides.fasta.genemark" -a -d
#FGS predict
echo "running frag gene scan"
#run_FragGeneScan.pl -genome=$out_dir"contigsWithPeptides.fasta" -out=$fgs_out"contigsWithPeptides.fasta.FGS" -complete=1 -train=complete -thread=$n_cores
#prodigal predict
echo "running prodigal"
#prodigal -i $out_dir"contigsWithPeptides.fasta" -o $prodigal_out"contigsWithPeptides.fasta.prodigal.ffn" -a $prodigal_out"contigsWithPeptides.fasta.prodigal.faa"
#post process genemarks output file and split between the prtein sequences and gene sequences
#python3 processGeneMark.py $genemark_out"contigsWithPeptides.fasta.genemark"

echo "constructing blast databases"
#construct blastDB
#makeblastdb -in $genemark_out"contigsWithPeptides.fasta.genemark.faa" -dbtype prot
#makeblastdb -in $fgs_out"contigsWithPeptides.fasta.FGS.faa" -dbtype prot
#makeblastdb -in $prodigal_out"contigsWithPeptides.fasta.prodigal.faa" -dbtype prot



if [ ! -d $blast_out ]; then
    mkdir $blast_out
fi

echo "blasting..."
#do protein blast
echo "blasting it against genemark database"
#blastp -db $genemark_out"contigsWithPeptides.fasta.genemark.faa" -query $out_dir"allExpressedORFsSoFar_c1.0.fasta" -num_threads $n_cores -outfmt 7 -out $blast_out"allExpressedORFsSoFar_c1.0_genemark_blastout.txt"
echo "blasting it against FragGeneScan database"
#blastp -db $fgs_out"contigsWithPeptides.fasta.FGS.faa" -query $out_dir"allExpressedORFsSoFar_c1.0.fasta" -num_threads $n_cores -outfmt 7 -out $blast_out"allExpressedORFsSoFar_c1.0_FGS_blastout.txt"
echo "blasting it against prodigal database"
#blastp -db $prodigal_out"contigsWithPeptides.fasta.prodigal.faa" -query $out_dir"allExpressedORFsSoFar_c1.0.fasta" -num_threads $n_cores -outfmt 7 -out $blast_out"allExpressedORFsSoFar_c1.0_prodigal_blastout.txt"

#process the blast outputs and create a summary table
python3 summarize_blastp.py  $blast_out"allExpressedORFsSoFar_c1.0_FGS_blastout.txt" $blast_out"allExpressedORFsSoFar_c1.0_prodigal_blastout.txt" $blast_out"allExpressedORFsSoFar_c1.0_genemark_blastout.txt" $expressedORFs_dir
#create a protein database using the refseq predicted proteomes
python3 get_refseq_proteomes.py $blast_out"blastHitsSummary.tsv" $contigs2genomes_f $proteomes_dir $out_dir
#create a protein blast database using the proteomes of the reference sequences
makeblastdb -in $out_dir"refseq_out/contigsWithPeptides.fasta.refseq.faa" -dbtype prot
#blast the ORFs with 100% hits against the refseq proteome database 
blastp -db $out_dir"refseq_out/contigsWithPeptides.fasta.refseq.faa" -query $out_dir"ORFsWithPeptidesWith100%hit.fasta" -num_threads $n_cores -outfmt 7 -out $blast_out"ORFsWithPeptidesWith100%hit_refseq_blastout.txt"
