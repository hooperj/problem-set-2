#! /usr/bin/env bash 

#These files are in the `data-sets/` directory.

#- Fasta file with human genome sequence: `fasta/hg19.chr22.fa`
#- BED file containing ChIP-seq peaks for H3K4me3 in Hela cells:
#  `bed/encode.h3k4me3.hela.chr22.bed.gz`
#- BED file with all genes in hg19: `bed/genes.hg19.bed.gz`.
#- File containing peak calls for ENCODE transcription factor ChIP-seq
#  experiements: `bed/encode.tfbs.chr22.bed.gz`.
#- Bedgraph with CTCF ChIP-seq data in bedGraph format: `bedgraph/ctcf.hela.chr22.bg.gz`
#- A "genome file" with chromosome size info: `genome/hg19.genome`
#- A file containing transcription start sites (TSS) for `chr22`: `bed/tss.hg19.chr22.bed.gz`


datadir=/WorkFiles/courses/JaysClass/data-sets

# create an array to store answers in and then use for loop to generate 
#	the properly formatted outputs. Don't forget that these arrays
#	are zero-indexed so that the first element is addressed with [0]

a=()
#1. identify the size of the largest overlap between CTCF and H3K4me3 locations.
a[1]=$(intersectBed -wao -a $datadir/bed/encode.h3k4me3.hela.chr22.bed.gz\
		  -b $datadir/bed/encode.tfbs.chr22.bed.gz |\
	sort -k 15 -g |\
# or sort -gk15  or -k15g but not -kg15
    	tail -n 1 |\
    	awk '{print $15}' )
# or cut -f 15)
#echo 'answer-1:' ${a[1]}


#2. Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#      19,000,500 on chr22 of `hg19` genome build. Report the GC content
#      as a fraction (e.g., 0.50).
#a[2]=$(gzcat $datadir/fasta/hg19.chr22.fa.gz|\
#	grep -v ^\>chr|\
#	tr '\n' ' '|sed 's/ //g' |\ 
#	cut -c 19000000-19000499|\
#	awk '{ gsub("[^C|G]", ""); print length/500 }' ) 
#echo 'answer-2: 0.384'
a[2]=0.384


#3. Use BEDtools to identify the length of the CTCF ChIP-seq peak in the
#	 bed/encode.tfbs.chr22.bed.gz file (i.e., interval) that has the 
#	largest mean signal in ctcf.hela.chr22.bg.gz.
a[3]=$(bedtools map -o mean -c 4 -b $datadir/bedtools/ctcf.hela.chr22.bg.gz\
	 -a $datadir/bed/encode.tfbs.chr22.bed.gz|\
	sort -k5nr| head -n 1|\
	awk 'BEGIN {FS="\t"} {print $3-$2}' )
#echo 'answer-3.1:' ${a[3]}

#4. Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of
#	a TSS) with the highest median signal in `ctcf.hela.chr22.bg.gz`.  Report
#	the gene name (e.g., 'ABC123')
# need to bedtools map -c 4 -o median the ctcf peaks in  -b bedtools/ctcf.hela.chr22.bg.gz
#	onto the  TSS in -a  bed/tss.hg19.chr22.bed.gz and then ...
a[4]=$(mapBed -c 4 -o median -a $datadir/bed/tss.hg19.chr22.bed.gz\
	 -b $datadir/bedtools/ctcf.hela.chr22.bg.gz|\
	sort -k7nr|\
	head -n1 |\
	cut -f 4)
#echo 'answer-4:' ${a[4]}


#5. Use BEDtools to identify the longest interval on `chr22` that is not
#	covered by `genes.hg19.bed.gz`. Report the interval like `chr1:100-500`.
a[5]=$(gzcat $datadir/bed/genes.hg19.sorted.bed.gz |\
	 bedtools complement -i - -g $datadir/genome/hg19.genome|\
	awk 'BEGIN {OFS="\t"} {if($1=="chr22") print $0, $3-$2}'|\
	sort -k4nr| head -n1|\
	awk 'BEGIN {OFS=""} { print $1, ":", $2, "-", $3}')
#echo 'answer-5:' ${a[5]} 

#6. Use one or more BEDtools that we haven't covered in class.
#  Take my 4992 'newgenes.bed'  converted to bed format from exon+transcript Stringmerge.gtf;
#	a gtf of known ensemble genes/transcripts/exons/etc.
#  	a bed file of 'repeats' combining repeatMasker, pseudogenes, that nasty LSU-rRNA_Hsa in WorkFiles/InfoMethodResource/mm10/DF0000772.mm10.hits.gz- to be built: 
#	a bed file of ESTs from UCSC table browser
# First use bedtools intersect (stranded) find exons that overlap known exons (minimum overlap 30 nt): then use awk to make a .bed output file of 'known' exons and transripts and of 'stillNovel.bed' exons and tanscripts
# then use betools subtract (not stranted) on 'stillNovel.bed' and 'junk.bed' to find unique EXON fragments- that do not overlap junk
# then use awk to filter for unique exon fragments longer than 30 nt and to write those exons and transcripts into 'novel.repeat.bed' and into 'STILLNovel.bed'. Note: aa transxript with ANY unique exon fragement longer than 30 nt is a novel transript
# then use  bedtoools interset -wao -s (stranded) on STILLNovel.bed and EST.bed to find transripts that overlap known ESTs (min overlap % ???) and write those exons and transripts to EST.novel
# same shit, but now take utterly novel - bedtools intersect -v to utterlyNovel.bed

 


end=(${!a[@]})   # put all the indices in an array
end=${end[@]: -1}    # get the last one
#echo 'length of answer: '$end

for i in 1 2 3 4 5 
    do
        echo 'answer-'$i':' ${a[$i]} 
    done
