# Detecting Human Pathogenic Bacteria GSM in Environmental Metagenomics

The characteristic sequences of human pathogenic bacteria are provided in FASTA format, along with an example workflow for alignment and statistics against environmental shotgun metagenomic sequences.

## Tools to be prepared:

Trimmomatic-0.39

perl

mmseqs

R

## 1. Quality Control

For pair-end sequencing, all samples should be placed in the same folder, and the sequencing results for each end should be named with _1 and _2. Perform quality control on all fq format files in the folder, and save the QC results with the original file name followed by the _p suffix in the ./qc folder.

Here is an example of the original sequence files for the samples:

>ls *.fq

>test0409A_1.fq  test0409A_2.fq  test0431A_1.fq  test0431A_2.fq

>mkdir qc

>for seq in `ls *.fq`;do
fname=`echo ${seq%%.fq}`
java -jar yourdir/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -phred33 -trimlog trim.log ${seq} ./qc/${fname}_p.fq LEADING:3 TRAILING:3 MINLEN:50 SLIDINGWINDOW:4:20
done



## 2. Calculate Normalized Values

Using Perl, calculate the correction values for the GSM matching counts from the quality-controlled sequences and save the file as stdnum.txt. The matching counts of GSM with the sequences need to be divided by the normalization value to correct and obtain N_st.

The formulas are as follows:

>N_st=(N_gsm×10^10)/((l_1-50)×n_1+(l_2-50)×n_2 )  (pair-end sequencing)

>N_st=(N_gsm×10^10)/((l-50)×n) （single-end sequencing）

N_st represents the normalized GSM counts. l_1 and l_2 denote the average read lengths of the paired-end sequencing data, while n_1 and n_2 represent the average read counts of the paired-end sequencing data. The value of N_st signifies the probability of detecting a selected pathogenic GSM at every 10^10 positions of 50 base pairs, each within the raw sequence dataset. 


### for pair-end sequencing:

>for file in `ls ./qc/*_p.fq`;do
smpname=`echo ${file} | sed 's|.*/\(.*\)_p\.fq|\1|'`
avglen=`perl -ne 'BEGIN{$min=1e10;$max=0;}next if ($.%4);chomp;$read_count++;$cur_length=length($_);$total_length+=$cur_length;END{print $total_length/$read_count,qq{ \n}}' ${file}`
linenum=`wc -l ${file}|awk 'BEGIN{FS=" "}{print $1}'`
readsnum=$((${linenum}/4))
echo -e "${smpname}\t${avglen}\t${readsnum}" >> seqstat.txt
done

>/usr/bin/Rscript /home/litong/2023dec/workflow/stdnum.R


### for single-end sequencing:

>echo -e "smpname\tstdnum" > stdnum.txt
for file in `ls ./qc/*_p.fq`;do
smpname=`echo ${file} | sed 's|.*/\(.*\)_p\.fq|\1|'`
avglen=`perl -ne 'BEGIN{$min=1e10;$max=0;}next if ($.%4);chomp;$read_count++;$cur_length=length($_);$total_length+=$cur_length;END{print $total_length/$read_count,qq{ \n}}' ${file}`
linenum=`wc -l ${file}|awk 'BEGIN{FS=" "}{print $1}'`
readsnum=$((${linenum}/4))
stdnum=$(echo "scale=4;(${avglen}-50)*${readsnum}/10^10" | bc)
echo -e "${smpname}\t${stdnum}" >> stdnum.txt
done


## 3. Alignment of GSM with Sample Sequences

Use the mmseqs to align pathogenic bacteria GSM with environmental sample sequences, considering only perfect matches. In the gsmfa folder, there are 10 fasta files, each containing characteristic sequences from 239 taxa of pathogenic bacteria, with 100 sequences for each species. When calculating, take the average of the alignment results from these 10 fasta files.

>mkdir results

>for i in `seq 1 10`;do
for file in `ls ./qc/*_p.fq`;do
outname=`echo ${file} | awk -F'/' '{sub("_p.fq", "", $NF); print $NF}'`
mmseqs easy-search ./gsmfa/gsmset${i}.fa ${file} ./results/${outname}_set${i}.out ./tmp --search-type 3 --min-aln-len 50 --min-seq-id 1 
done
done


Results are stored in ./results.

## 4. Result Counting and Normalization

Count the matching numbers between GSM and samples, and correct them according to the previously calculated normalization values.

--rdir: Folder where matching results are located.

--thr: If the number of detected GSM species within the same species is less than this value, it is not included in the detection results. Recommended value is 2.

--speid: The numerical prefix before the GSM identifier represents the species name (this file has been provided).

--outdir: File where the output results are stored.

>/usr/bin/Rscript stat_and_std.R --rdir ./results --thr 2 --speid ./speid.txt --outdir ./out --std ./stdnum.txt


Below is an example of a result file:

>cat ./out/PathogensStat.txt

>"spename"	"test0409A"	"test0431A"
"Clostridioides difficile"	32.2751839816847	13.8916324234751
"Clostridium novyi"	64.5503679633694	29.0755097235525
"Eggerthella lenta"	1.29100735926739	5.81510194471049
"Listeria monocytogenes"	0	2.26142853405408

The first column is the name of the pathogen species, the column name is the name of each sample, and the value is the normalized N_st.
