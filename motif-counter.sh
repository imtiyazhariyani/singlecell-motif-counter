#!/bin/bash
#SBATCH --job-name=FCA_SdicCent_MotifCounter           ## Name of the job. !!CHANGE THIS TO APPROPRIATE NAME!!
#SBATCH -p standard                                     ## partition/queue name
#SBATCH --nodes=1                                       ## use 1 node, don’t ask for multiple
#SBATCH --ntasks=4                                      ##number of CPUS on 1 node
#SBATCH --time=72:00:00                                 ##time limit for job
#SBATCH --mem-per-cpu=8G                                ##memory put CPU; if ntasks = 4, total mem = 32Gb
#SBATCH --output=htstream_%A_%a.out         ## File to which STDOUT will be written ## slurm error file, %x - job name, %A job id
#SBATCH --error=htstream_%A_%a.err          ## File to which STDERR will be written ## slurm output file, %x - job name, %A job id

set -e
set -u
set -o pipefail

##This script was repeated for both L3 and FCA adult datasets for Sdic-aggregate, Sdic paralogs and sw using known core and extended motifs
##Below is an example for Sdic_1 on the FCA dataset

##Inpath and outpath, should take processed reads from 01- and put the outputs into 03- - CHANGE FOR EACH TISSUE

INPATH="/share/crsp/lab/jranz/ihariyan/motif-counter/FCA/HTS_PreprocDUP"
OUTPATH="/share/crsp/lab/jranz/ihariyan/motif-counter/FCA/Centered"

##R1 Variables (Rev Compl) - CHANGE FOR EACH GENE
motifID_cur="Sdic_A"     ## current motif ID, eg. Sdic1. - CHANGE FOR EACH GENE (SAME FOR R1&R2)
coremotifSeq_cur1="CCCATCACCACGCTATTGAT"    ## the seq of current core motif. - CHANGE FOR EACH GENE    
refSeq_cur1="GACCCCGGAGCGCAGGCCGTGGCGCGAGGCGGAGTAGACGTAGCCGTCCTCACTGCCCATCACCACGCTATTGATCTCGTTGGCCGGGAAGGCCATCGATGTAATGGCAATGGCCTTAGACTGGCGCTGC" ##!!the seq of current extended motif. - CHANGE FOR EACH GENE
LCore_cur1="20" ## !!the length of current core motif seq - CHANGE FOR EACH GENE (SAME FOR R1&R2)
LRef_cur1="130" ## the length of current motif seq. - keep constant at 130 (same for R1&R2)
SRef_cur1="56"  ## !!the start pos of current core motif in ref. - CHANGE FOR EACH GENE
ERef_cur1="75"  ## !!the end pos of current core motif in ref. - CHANGE FOR EACH GENE

##R2 Variables - CHANGE FOR EACH GENE
coremotifSeq_cur2="ATCAATAGCGTGGTGATGGG"    ## the seq of current core motif. - CHANGE FOR EACH GENE
refSeq_cur2="GCAGCGCCAGTCTAAGGCCATTGCCATTACATCGATGGCCTTCCCGGCCAACGAGATCAATAGCGTGGTGATGGGCAGTGAGGACGGCTACGTCTACTCCGCCTCGCGCCACGGCCTGCGCTCCGGGGTC" ## !!the seq of current extended motif. - CHANGE FOR EACH GENE
LCore_cur2="20" ## !!the length of current core motif seq - CHANGE FOR EACH GENE (SAME FOR R1&R2)
LRef_cur2="130" ## the length of current motif seq. - keep constant at 130
SRef_cur2="56"  ## !!the start pos of current core motif in ref. - CHANGE FOR EACH GENE
ERef_cur2="75"  ## !!the end pos of current core motif in ref. - CHANGE FOR EACH GENE

##DO NOT MODIFY PAST HERE

for file in "${INPATH}"/*; do
        SAMPLE=$(basename "$file" | sed 's/_R2.*//')
#R1
        zcat ${INPATH}/${SAMPLE}_R2_001.fastq.gz_SE.fastq.gz | grep -A 2 -B 1 ${coremotifSeq_cur1} - | sed '/^--/d' | gzip - > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1.fastq.gz
        zcat ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1.fastq.gz | grep -A 1 '^@'     | sed '/^--/d'     | sed 's/ /_/g'     | awk '/^@/&&NR>1{print "";}{ printf "%s",/^@/ ? $0"\t":$0 }' | awk -v n=${LCore_cur1} -v b=${coremotifSeq_cur1}  -v L=${LRef_cur1} -v S=${SRef_cur1} -v E=${ERef_cur1} -v RSeq=${refSeq_cur1} -v LRead=98 '{SRead=index($2,b);ERead=SRead+n-1;HeadLength=SRead-1;if(HeadLength>S-1)HeadLength=S-1; TailLength=LRead-ERead;if(TailLength>L-E)TailLength=L-E; SRefOverlap=S-HeadLength; ERefOverlap=E+TailLength; SReadOverlap=SRead-HeadLength; EReadOverlap=ERead+TailLength; LOverlap=EReadOverlap-SReadOverlap+1; SeqRefOverlap=substr(RSeq,SRefOverlap,LOverlap) ; SeqReadOverlap=substr($2,SReadOverlap,LOverlap) ; print $0"\t"LRead"\t"SRead"\t"ERead"\t"HeadLength"\t"TailLength"\t"SRefOverlap"\t"ERefOverlap"\t"SReadOverlap"\t"EReadOverlap"\t"LOverlap"\t"SeqRefOverlap"\t"SeqReadOverlap}' > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1_overlap
        awk '{if($3!=-1 && $12>=70)print}' ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1_overlap > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1_Reads
        awk '{if($3!=-1 && $12>=70 && $13==$14)print}' ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1_overlap > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R1.PerfReads

#R2
        zcat ${INPATH}/${SAMPLE}_R2_001.fastq.gz_SE.fastq.gz | grep -A 2 -B 1 ${coremotifSeq_cur2} - | sed '/^--/d' | gzip - > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2.fastq.gz
        zcat ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2.fastq.gz | grep -A 1 '^@'     | sed '/^--/d'     | sed 's/ /_/g'     | awk '/^@/&&NR>1{print "";}{ printf "%s",/^@/ ? $0"\t":$0 }' | awk -v n=${LCore_cur2} -v b=${coremotifSeq_cur2}  -v L=${LRef_cur2} -v S=${SRef_cur2} -v E=${ERef_cur2} -v RSeq=${refSeq_cur2} -v LRead=98 '{SRead=index($2,b);ERead=SRead+n-1;HeadLength=SRead-1;if(HeadLength>S-1)HeadLength=S-1; TailLength=LRead-ERead;if(TailLength>L-E)TailLength=L-E; SRefOverlap=S-HeadLength; ERefOverlap=E+TailLength; SReadOverlap=SRead-HeadLength; EReadOverlap=ERead+TailLength; LOverlap=EReadOverlap-SReadOverlap+1; SeqRefOverlap=substr(RSeq,SRefOverlap,LOverlap) ; SeqReadOverlap=substr($2,SReadOverlap,LOverlap) ; print $0"\t"LRead"\t"SRead"\t"ERead"\t"HeadLength"\t"TailLength"\t"SRefOverlap"\t"ERefOverlap"\t"SReadOverlap"\t"EReadOverlap"\t"LOverlap"\t"SeqRefOverlap"\t"SeqReadOverlap}' > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2_overlap
        awk '{if($3!=-1 && $12>=70)print}' ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2_overlap > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2_Reads
        awk '{if($3!=-1 && $12>=70 && $13==$14)print}' ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2_overlap > ${OUTPATH}/${motifID_cur}.${SAMPLE}.R2.PerfReads


done

#Annotate the identfied perfect matches above as reads belonging to a previously identified cell type

module load anaconda
source activate scTE

samples=("FCA59_Male_testis_adult_1dWT_Fuller_sample1_S44_L003" \
         "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S48_L003" \
         "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S49_L003" \
         "FCA59_Male_testis_adult_1dWT_Fuller_sample1_S56_L003" \
	 	 "FCA60_Male_testis_adult_1dWT_Fuller_sample2_S59_L004" \
	 	 "FCA60_Male_testis_adult_1dWT_Fuller_sample2_S61_L004" \
	 	 "FCA60_Male_testis_adult_1dWT_Fuller_sample2_S67_L004" \
	 	 "FCA60_Male_testis_adult_1dWT_Fuller_sample2_S70_L004" \
	 	 "FCA61_Male_testis_adult_1dWT_Fuller_sample3_S60_L004" \
	 	 "FCA61_Male_testis_adult_1dWT_Fuller_sample3_S65_L004" \
	 	 "FCA61_Male_testis_adult_1dWT_Fuller_sample3_S68_L004" \
         "FCA61_Male_testis_adult_1dWT_Fuller_sample3_S69_L004")

for sample in "${samples[@]}"; do
	python v15_cellID_annotate.py \
		--input-UMI-fastq "/pub/ihariyan/ranz_lab/data/single-cell/fly-cell-atlas/${sample}_R1_001.fastq.gz" \
    	--input-perf-reads "/share/crsp/lab/jranz/ihariyan/motif-counter/FCA/Centered/Sdic_1.processed.${sample}.R2.PerfReads" \
    	--output "${sample}.Sdic_1.R2.out"
done

#python script is provided separately

cat *.Sdic_1.R1.out > Sdic_1.R1.out

cat *.Sdic_1.R2.out > Sdic_1.R2.out

cut -f2 Sdic_1.R1.out | sort | uniq -c > Sdic_1.R1.counts

cut -f2 Sdic_1.R2.out | sort | uniq -c > Sdic_1.R2.counts

