#!/bin/bash
echo -e "
    __  __ ______ ____   ______ __  __ __     ______ _____
   / / / // ____// __ \ / ____// / / // /    / ____// ___/
  / /_/ // __/  / /_/ // /    / / / // /    / __/   \__ \ 
 / __  // /___ / _, _// /___ / /_/ // /___ / /___  ___/ / 
/_/ /_//_____//_/ |_| \____/ \____//_____//_____/ /____/  
                                                          
"

echo ""
echo "High-throughput Epidemiological Reconstruction and Clustering for Uncovering Lineages from Environmental SARS-CoV-2."
echo ""
echo "Quality:"${1}"/Noise Cutoff:"${2}
echo "Analyzing Spike gene from nt "${3}" to "${4}
echo "Read size between "${5}" and "${6} "nt"
echo "Trimming "${8}"nt"
echo "Additional sites "${9}
echo ""
echo "Visit us at https://github.com/garcia-nacho/HERCULES"
sleep 5s
echo ""
echo ""

RefBowtie2=/home/docker/CommonFiles/reference/SpikeSars-CoV-2
RefSpike=/home/docker/CommonFiles/reference/SpikeRef.fa
Tools=/home/docker/CommonFiles/Tools

basedir=/Data

mkdir -p /home/docker/results/
Rscript /home/docker/CommonFiles/Tools/FolderRenamer.R

source activate nextclade
echo "Downloading Nextclade database"
nextclade dataset get --name 'sars-cov-2' --output-dir '/home/docker/nc_sars-cov-2'
conda deactivate

#Legacy
#echo ""
#echo "Preparing Spike References"
#cat /home/docker/CommonFiles/Variants/*.fa* > /Data/variantRefs.fasta
#Rscript /home/docker/CommonFiles/Tools/Spike_Extractor.R
#rm /Data/variantRefs.fasta

for dir in $(ls -d */)
do
    
    numberoffiles=$(ls ${dir}*.fastq.gz | wc -l)
    SKIP="FALSE"

    if (( ${numberoffiles} == 1 ))
    then
    FILESIZE=$(stat -c%s ${dir}/*.fastq.gz)

    if (( ${FILESIZE} < 1000000))
    then
    SKIP="TRUE"
    echo "Skipping "${dir}
    fi
    fi

    if [[ ${numberoffiles} > 1  || ${SKIP} == "FALSE" ]] 
    then
    echo ""
    echo "Processing "${numberoffiles}" fastq.gz files in "${dir}
    echo ""
    
    

    cd ${dir}
	  Reads=$(ls *.fastq.gz)
    cat ${Reads} > ${dir%/}.fastq.gz
    gzip -d ${dir%/}.fastq.gz
    seqkit seq ${dir%/}.fastq -M ${6} -m ${5} -Q ${1} > ${dir%/}.filtered.fastq

    if (( ${8} != 0))
    then
    echo "Trimming "${8}"nt"
    source activate cutadaptenv
    cutadapt --cut ${8} -o ./trimmed1.fastq ./${dir%/}.filtered.fastq 
    rm ./${dir%/}.filtered.fastq
    cutadapt --cut -${8} -o ./trimmed2.fastq ./trimmed1.fastq
    rm ./trimmed1.fastq
    mv ./trimmed2.fastq ./${dir%/}.filtered.fastq
    conda deactivate
    fi

    rm ${dir%/}.fastq
    #(bowtie2 -p 8 -x ${RefBowtie2} -U ${dir%/}.fastq -S ${dir%/}.sam) 2> ${dir%/}_Bowtie2summary.txt
    #(tanoti -r ${RefSpike} -i ${dir%/}.fastq -o ${dir%/}.sam -u) 2> Tanoti_${dir%/}.summary.txt
    minimap2 --secondary=no -ax map-ont ${RefSpike} ${dir%/}.filtered.fastq > ${dir%/}.sam  
    samtools view -F 1024 -F 256 -F4 -F 2048 -bS ${dir%/}.sam | samtools sort -o ${dir%/}.sorted.bam
    samtools index ${dir%/}.sorted.bam
    
    #Legacy
    #minimap2 --secondary=no -ax map-ont /Data/VariantSpike.fasta ${dir%/}.filtered.fastq > ${dir%/}_voc.sam  
    #Rscript /home/docker/CommonFiles/Tools/SamParser.R

    ls
    samtools mpileup -aa -A -d 0 -q 0 --reference ${RefSpike} ./${dir%/}.sorted.bam | ivar consensus -t 0 -n N -m 20 -p ${dir%/}_consensus

    rm ${dir%/}.sam
    #Legacy
    #rm ${dir%/}_voc.sam
    cat ${dir%/}.filtered.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ${dir%/}_read_length.txt

    #rm ${dir%/}_voc.sorted.bam
    #rm ${dir%/}_voc.sorted.bam.bai

    #Legacy
    #echo "Running AF search"
    #echo ""
    #cp /home/docker/CommonFiles/reference/reference.msh ./reference.msh
    #cp /home/docker/CommonFiles/Tools/mash ./mash
    
    #seqkit fq2fa ${dir%/}.filtered.fastq -o ${dir%/}.uncompressed.fasta
    rm ${dir%/}.filtered.fastq
    
    #Rscript /home/docker/CommonFiles/Tools/MashRunner_V2.R
    mv ${dir%/}_read_length.txt /home/docker/results/${dir%/}_read_length.txt
    #mv ${dir%/}_AF.lineages.csv /home/docker/results/${dir%/}_AF.lineages.csv
    #rm ${dir%/}.uncompressed.fasta
    #rm ./mash
    #rm ./reference.msh

    ${Tools}/FINex2 -f ${dir%/}.sorted.bam > ${dir%/}.noise.tsv 
    cat ${dir%/}_consensus.fa ${RefSpike} > spike.cons.fa
    
    source activate nextclade
    nextclade run --input-dataset /home/docker/nc_sars-cov-2 --output-csv dummy.csv --output-fasta spike.cons.aligned.fa spike.cons.fa
    conda deactivate
    #nextalign -i spike.cons.fa -o spike.cons.aligned.fa -r /home/docker/CommonFiles/reference/SpikeRef.fa

    Rscript ${Tools}/AnalysisWW_V03.R $(pwd)/ ${2} ${3} ${4}

    mv Coverage.pdf /home/docker/results/${dir%/}.coverage.pdf
    mv Noise.pdf /home/docker/results/${dir%/}.noise.pdf
    mv ${dir%/}.noise.tsv  /home/docker/results/${dir%/}.noise.tsv 
    mv ${dir%/}.sorted.bam /home/docker/results/${dir%/}.sorted.bam
    mv ${dir%/}.sorted.bam.bai /home/docker/results/${dir%/}.sorted.bam.bai
    mv *_consensus.qual.txt /home/docker/results/${dir%/}_consensus.qual.txt 
    mv ${dir%/}_consensus.fa /home/docker/results/${dir%/}_consensus.fa
    mv spike.cons.aligned.fa /home/docker/results/${dir%/}_spike.cons.aligned.fa
    #mv ${dir%/}_voc_depth.csv /home/docker/results/${dir%/}_voc_depth.csv    

    rm dummy.csv
    rm Rplots.pdf
    rm spike.cons*
    rm 
    cd ${basedir}  
    fi

done


cp -R /home/docker/results /Data/results

    if [ -d "/Previous" ]
    then
    echo "======================="
    echo "Integrating old results"
    echo "======================="
    mkdir /Data/results/OldResults
    cp /Previous/* /Data/results/OldResults
    fi

cd /Data/results

mkdir bam analysis QC sequences
mkdir analysis/SankeyPlots

mv /Data/results/*.bam /Data/results/bam
mv /Data/results/*.bai /Data/results/bam
mv /Data/results/*.fa /Data/results/sequences
mv /Data/results/*.qual.txt /Data/results/QC
mv /Data/results/*noise.tsv /Data/results/QC
mv *_consensus.qual.txt /Data/results/QC 
mv *_read_length.txt /Data/results/QC 
mv *pdf /Data/results/QC

#Kmer search



Rscript /home/docker/CommonFiles/Tools/ParallelBamExtraction.R ${2} ${3} ${4} ${9} ${11}

mv *.html analysis
mv *Sankeyplot.Mutations.pdf analysis/SankeyPlots

#Legacy
#Rscript /home/docker/CommonFiles/Tools/MappedVariantPlotter.R
#Rscript /home/docker/CommonFiles/Tools/MashPlotter.R

#mkdir /Data/results/analysis/Legacy
#mv CountVariant_groupped_byVariant.pdf /Data/results/analysis/Legacy/CountVariantMapped_groupped_byVariant.pdf
#mv CountVariant_groupSample.pdf /Data/results/analysis/Legacy/CountVariantMapped_groupSample.pdf
#mv *tMash*.pdf /Data/results/analysis/Legacy
#mv VariantMappedReads.xlsx /Data/results/analysis/Legacy
#mv AF_lineages.xlsx /Data/results/analysis/Legacy/VariantMashReads.xlsx
#mv RatioVariant_groupSample.pdf /Data/results/analysis/Legacy/RatioVariantMapped_groupSample.pdf
#mv RatioVariant_groupVariant.pdf /Data/results/analysis/Legacy/RatioVariantMapped_groupVariant.pdf
#mv RatioVariantMapped_Stacked.pdf /Data/results/analysis/Legacy/RatioVariantMapped_Stacked.pdf

#mkdir /Data/results/analysis/Widgets 
mv /Data/results/analysis/*.html /Data/results/analysis/Widgets
cp /home/docker/CommonFiles/Dashboard.html /Data/results/analysis/Dashboard.html

mkdir /Data/results/analysis/SingleMutations
mv *SingleMutation* /Data/results/analysis/SingleMutations

mv /Data/VariantSpike.fasta /Data/results/analysis/ReferencesSpike.fasta
#mv /Data/results/*_voc_depth.csv /Data/results/QC


mv /Data/results/analysis/Clade_Barplot.pdf /Data/results/QC
#mv /Data/results/analysis/CountVariant_groupped_byVariant.pdf /Data/results/analysis/CountVariantMapped_byVariant.pdf
#mv /Data/results/analysis/CountVariant_groupSample.pdf /Data/results/analysis/CountVariantMapped_bySample.pdf 
mv /Data/results/analysis/CoverageMosaic.pdf /Data/results/QC
mv /Data/results/analysis/NoiseMosaic.pdf /Data/results/QC
mv /Data/results/analysis/Pangolin_Barplot.pdf /Data/results/QC
#mv /Data/results/analysis/RatioVariant_groupSample.pdf /Data/results/analysis/RatioVariantMapped_bySample.pdf
#mv /Data/results/analysis/RatioVariant_groupVariant.pdf /Data/results/analysis/RatioVariantMapped_byVariant.pdf  
mv /Data/results/analysis/ReferencesSpike.fasta /Data/results/QC
mv /Data/results/analysis/Variants.nextclade /Data/results/QC
#mv /Data/results/*AF.lineages.csv /Data/results/QC
mv /Data/results/*.gz /Data/results/QC
mv /Data/results/*read_length.txt /Data/results/QC
mv /Data/results/*.pdf /Data/results/analysis
mv /Data/results/*.xlsx /Data/results/analysis
rm /Data/results/analysis/Rplots.pdf
rm -rf /Data/results/OldResults
rm /Data/results/MSAFastas

mkdir /Data/results/WWDB
mv /Data/results/MSAFastas/* /Data/results/WWDB
rm -rf mv /Data/results/MSAFastas/
cp /Data/results/QC/*noise.tsv /Data/results/WWDB

if [[ ${10} != "0" ]]
then

cd /Data/results
Rscript /home/docker/CommonFiles/Tools/KmerSearchDockerParallel.R ${10} 
rm -rf kmeruncompressed
mv ExtractedKmer analysis/KmerSearch/
mkdir /Data/results/sequences/KmerSeqs
mv ExtractedKmer analysis/KmerSearch/
mv analysis/KmerSearch/* sequences/KmerSeqs

mv KmerSearchBarplot.pdf analysis/KmerSearch/KmerSearchBarplot.pdf  
mv KmerSearchResults.xlsx analysis/KmerSearch/KmerSearchResults.xlsx  

fi
