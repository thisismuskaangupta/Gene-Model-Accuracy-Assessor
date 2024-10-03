# Gene-Model-Accuracy-Assessor
this contains the code written for the MSc Bioinformatics Project. Grade: A3

# Project Summary 
Gene structural annotation is the process of finding and storing information on ‘features’ within a genetic sequence, including but not limited to exons, introns, and promoters. These vary widely in quality and accuracy due to many factors, including how well-studied the organism is. Naturally, quality is lower for neglected divergent organisms. VEuPathDB, a bioinformatics resource centre, hosts information on about 700 neglected organisms, and there is a need to highlight to users when an available canonical model may be incorrect. There exists current software to judge the accuracy of gene models, such as MAKER2 and TRaCE, however these require computationally-expensive storage and preprocessing of RNA-seq read files, including trimming, alignments, and assembly. VEuPathDB outsources this processing to the European Bioinformatics Institute (EBI), and doesn’t store these files. In the manifest of files returned by EBI and stored on VEuPathDB servers is a junction dataset, derived from splice-split read alignments, containing information on location of junctions and the number of reads split across said junction. The primary objective of this program was to mirror the functionality of a metric called Annotation Edit Distance (AED) used in MAKER2 and TRaCE, using the existing junction data in lieu of computationally expensive transcript assemblies. Further, this project aimed to allot flags or comments for each gene to elucidate the validity of reference gene models supplied, and in the case that cumulative data from experiments suggests the existence of other models, to provide these and estimate a ‘best alternative’ transcript. 

A tool was developed that would use the junction dataset, a set of annotations, and a genome sequence file to analyse and comment on each gene model (transcript) in the annotation, and categorize this under a ‘case’, which would indicate additional information such as if the gene was exon-only, if flanking introns were found, if the transcript’s translation was valid, if it was supported by evidence from the junction dataset, etc. Further, gene models were predicted from junction evidence and compared to the canonical transcript(s) using AED, and a ‘best alternative’ was suggested. 

The tool, developed using _Toxoplasma gondii_ ME49 strain annotations b65 (old) and b68 (new), was tested on _Mucor lusitanicus_ CBS 277.49 strain annotations b60 (old) and b66 (new). It was able to capture improvements between older and updated genome annotations; it also provided a broad overlook on the differences between overall quality in the annotation of the two organisms – most genes in the updated annotations for _T. gondii_ were likely valid, while for _M. lusitanicus_, most were unsure. This was in line with the fact that _Mucor_ is much less studied, with lesser sequencing data available. A few different types of genes in _T. gondii_ were manually selected and their output looked at for a more microscopic evaluation. 

Several limitations were discussed and framed in the context of scope for future improvement.

# Command Line 
\< usage: summer project [-h] --junctions JUNCTIONS --gff_whole GFF_WHOLE --gff_introns_only GFF_INTRONS_ONLY --genome_sequence GENOME_SEQUENCE [--AED_threshold AED_THRESHOLD] [--intron_cutoff INTRON_CUTOFF] \>  

# Inputs 
For an organism of choice on VEuPathDB: 
  1. A folder containing all the junction data files across the whole genome. They may be organized in any way. 
  2. A GFF-format annotation file that contains features across the whole genome (including ‘intron’ features). Some newer GFF formats don’t contain intron features, they may conveniently be added using NBISweden’s AGAT toolkit: agat_sp_add_introns.pl tool, version: v1.0.0. 
  3. Another GFF-format file containing only intron features across the whole genome. This may be easily created from the previous file using grep on a Linux command line. (Following is an example command for the same.)
     
  \< grep intron whole_genome_with_introns.gff \> introns_only.gff \> 
  
  4. A fasta-format file containing the genome sequence. 
  5. (Optional) A threshold for AED, beyond which gene model predictions won’t be recorded. Default: 0.5 
  6. (Optional) A unique read depth cutoff for noisy junction data (with respect to the intron with the highest read depth for its respective gene). Default: 0.2. Junctions that have a read depth of less than (the cutoff * maximum read depth) will be removed. This rule is followed unless the read depth of the intron with the highest read depth is more than 5000. If it is more than 5000, then any intron with a read count of 1000 or more is selected.
     
# Outputs 
  1. A log file called ‘log.log’ containing information on –  
    a. total junctions in junction files, how many were found in annotations, how many could not be found in annotations, and how many introns were in annotations but not in junctions files;  
    b. reporting the creation of each output file once it is created;  
    c. reporting the number of genes found across three major categories of flags – ‘likely valid’, ‘likely invalid’, and ‘unsure’. 
  2. A GFF-format file called ‘all_introns.gff’ containing information on all introns in the junction files. 
  3. GFF-formatted annotations of predicted gene models appended to the supplied whole genome GFF file. 
  4. Two histograms called 'readsonly_depth.png' and 'readsandanns_depth.png' which demonstrate the distribution of unique read depth for introns found in junction files only (not annotated) and introns that were found in annotations respectively. 
  5. A violin plot called ‘cutoff.png’ demonstrating the distribution of the ratio of each intron’s read depth to the read depth of the intron with maximum read depth for the gene that it lies in, across the whole genome. 
  6. A file called 'gene_info.csv’ which is an information catalogue containing gene model accuracy data about all the genes. It contains the following columns:  
    a. Gene  
    b. Flag  
    c. Case 
    d. Ref Transcript  
    e. AED Scores (a semicolon-delimited list in the format - predictionID : AED;)
    f. Total Introns Found  
    g. Number of Valid Introns 
    h. Mean Read Depth Before Filtering  
    i. Mean Read Depth After Filtering  
    j. Maximum Read Depth  
    k. Cutoff Used
    l. Evidence of Validity of Ref Transcript?  
    m. Evidence of Existence of Alternative Transcript?  
    n. Best Alternative Model 
