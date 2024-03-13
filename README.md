# Wastewater SARS-CoV-2 Surveillance Pipeline
[![Docker](https://badgen.net/badge/icon/docker?icon=docker&label)](https://https://docker.com/)


SARS-CoV-2 wastewater surveillance using Nanopore long reads. Version 1.1

## Introduction
This pipeline aims to identify and quantify SARS-CoV-2 lineages at the read-level.   
<img src="/Overview.png" width="700">     

## Installation   
The recommended way to install and use the pipeline is by using the pre-built docker image
<code>docker pull ghcr.io/garcia-nacho/wastewater:v1.1</code>
However, you could download and edit this repo to create a custom analysis yourself: 
<code>git clone https://github.com/garcia-nacho/Wastewater_SARS-CoV-2/ </code>  
<code>docker build -t wastewater Wastewater_SARS-CoV-2</code>

## Run   
The pipeline accepts several parameters that are parsed to the image via the flag <code>-e parameter=value</code>

**Quality** <code>qual=value</code>  
The pipeline will filter the reads based on a quality cut-off.    
The default quality threshold is 15. But a quality limit of 10 also gives good results.    
A value of -1 will prevent any filtering.    

**Minimum Size** <code>m=value</code>  
The pipeline will filter reads with a size lower than the minimum value.
The default minimum size is 500nt

**Maximum Size** <code>M=value</code>
The pipeline will filter reads with a size higher than the maximum value.
The default minimum size is 1300nt

**Noise** <code>noise=value</code>
The pipeline will analyze the positions with a noise higher than the value provided by this parameter.
The default noise level is 0.15 but values between 0.05 and 0.10 could give good results.   
A very high noise parameter will cause fewer sites analyzed and therefore some mutations might be ignored.    
A very low noise parameter will cause the inclusion of more sites, which would cause a higher running time and might reduce the performance of the analysis.

**Start** <code>start=value</code>
This parameter defines where the sequenced region into the Spike gene starts.
The default start position is 1250   

**End** <code>end=value</code>
This parameter defines where the sequenced region into the Spike gene ends.
The default start position is 2250   

**Trimming** <code>trim=value</code>    
This parameter specifies if trimming is needed to remove the primers and how many nucleotides should be trimmed.
The default value is no trimming (trim=0)

**Positions of Interest** <code>poi="pos1,pos2,pos3"</code>   
This parameter defines which positions should be analyzed disregarding the noise level on them. The accept the standard mutation format considering the entire SARS-CoV-2 genome (e.g. "G22895C,T22896A,G22898A,A22910G')
The default value is empty.  

**Kmer**  <code>kmer="kmer1,kmer2,kmer3</code>   
This parameter allows you to screen for kmers in the reference and the reads. 
The default value is empty

The pipeline expects this folder structure:   
<pre>
./ExpXX         
  |-Sample1     
      |-File_XXXXX_1.fastq.gz       
      |-File_XXXXX_2.fastq.gz
      |-File_XXXXX_3.fastq.gz
      |-...
  |-Sample2      
      |-File_XXXXX_1.fastq.gz       
      |-File_XXXXX_2.fastq.gz
      |-File_XXXXX_3.fastq.gz
      |-... 
  |-...   

</pre>

The filename of the *.fastq.gz* files are irrelevant and the samples are named using the folder that containes them as name    

Note that older versions of docker might require the flag <code>--privileged</code> to run properly. 

## Integration of old results   
It is possible to integrate previous results of the pipeline to be analyzed again using the information of the new samples. 
To integrate the results of old analyses, a directory containing the files with the suffix *_b2f.tsv.gz* should be mapped into the /Previous folder of the docker container by using the following flag. 
<code>-v Path/to/database:/Previous</code>
   
   
## Running Examples   
Basic run using default settings:   
<code>docker run -it --rm -v $(pwd):/Data ghcr.io/garcia-nacho/wastewater:v1.1 </code>  
   
Read filter by a quality of 10:  
<code>docker run -it --rm -e qual=10 -v $(pwd):/Data ghcr.io/garcia-nacho/wastewater:v1.1 </code>   
   
Noise cut-off change to 0.1:       
<code>docker run -it --rm -e noise=0.1 -v $(pwd):/Data ghcr.io/garcia-nacho/wastewater:v1.1 </code>

Change the read size to allow reads between 100 and 2000nt    
<code>docker run -it --rm -e m=100 -e M=2000 -v $(pwd):/Data wastewater </code>   

Look into mutations of the BA.2.86 variant and use the previous results included in the directory WWDB   
<code>docker run -it --rm -v WWDB:/Previous -v $(pwd):/Data \ 
-e poi='G22895C,T22896A,G22898A,A22910G,C22916T,G23012A,C23013A,T23018C,T23019C,C23271T,C23423T,A23604G' \ 
ghcr.io/garcia-nacho/wastewater:v1.1</code>

Look into a kmer specific from the BA.2.86 variant and use the previous results included in the directory WWDB   
<code>docker run -it --rm -v WWDB:/Previous -v $(pwd):/Data -e kmer='TAAGCATAGTG' ghcr.io/garcia-nacho/wastewater:v1.1</code>


## Output   
The pipeline generates several folders: analysis, bam, QC, sequences and WWDB   
   
**analysis**    
Here you will find:   
-Sankey plots showing the most abundant combinations of mutations at nucleotide and amino acid levels   
-Barplots showing the relative abundance of the different combinations. Different granularity levels are included (amino acid sequence level, pangolin level, mutation-based-lineage)   
-Excel files with the raw results    
-A barplot showing the relative abundance of all the single mutations found in the samples.
-html widgets showing the distribution of the different lineages found in the different samples and in the entire dataset.

**bam**   
Here you will find the *bam* files containing the reads aligned against the Spike gene of the *wuhan-hu-1* strain   
   
**QC**   
Here you will find:   
-Plots showing the coverage of the samples over the spike gene   
-Noise at the different positions (Noise is defined here as ratio of bases not included in the consensus)   
-A noise.tsv file that contains the raw data regarding coverage and noise   
-consensus_qual.txt. Quality of the different bases called in the consensus   
   
**sequences**   
Here you find the fasta sequences of the consensus and variants found in each sample    

**WWDB**
Here you will find the files required to integrate the results into future analyses.
     
## Under the hood
The pipeline filters the reads according to a quality cut-off using *[seqkit](https://bioinf.shenwei.me/seqkit/)*. Reads are mapped against the reference Spike using *[minimap2](https://github.com/lh3/minimap2)* and filtered and sorted using *[samtools](http://www.htslib.org/)*. The resulting *bam* file is indexed using *[samtools](http://www.htslib.org/)*. The positions showing a mix of bases are identified using *[noisefinder](https://github.com/garcia-nacho/NoisExtractor)* and the bases connected with their read-ids are retreived using *[bbasereader](https://github.com/garcia-nacho/bbasereader)*. All positions are merged using the read-id as pivot column and the different variants are identified and saved in a fasta file. 
The pangolin lineages are assigned using the mutations on the references and reads. If there is no discrepancy between the mutations found in a read and the mutations found in a lineage, the read is assigned to that lineage. If there are discrepancies the reads are assigned to the closer lineage and the discrepancies are displayed. If there are several lineages close enough to the read in terms of mutations, all of them are displayed (Note that not all mutations are analyzed, only those dynamically decided and those defined with the *poi* flag).
The plots and analyses are generated in *[R](https://www.r-project.org/)* 
