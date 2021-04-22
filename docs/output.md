# Outputs

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes the data using the steps presented in the main README file.  
Briefly, its goal is to process single-cell ChIP-seq data from inDrop protocol.

The directories listed below will be created in the output directory after the pipeline has finished. 

### Barcode Matching

![MultiQC](images/barcode.png)

Barcodes are composed of three indexes originating from three different libraries. 

![MultiQC](images/scChIPseq_barcode_plot-1.png)

If 40% or more sequences are correctly barcoded ("Barcoded" in the legend), barcoding step went well. If it is between 20% and 40% be carrefull and if it less, barcoding step went wrong.

### Bowtie2

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner) software is used to aligned barcodes to index reference libraries. 

- **SE mapped uniquely** : successfully aligned sequence to a unique index.  
- **SE multimapped** : sequence aligned on multiple indexes.  
- **SE not aligned** : sequence corresponding to none index. 

A drop in successfull alignement percentage is observed after the first index. It is due to the sequential ligation of indexes (first the index 1 is ligated, then the 2, then the 3) that leads to a loss at each round.
It is attempted that arround 75% of the first index map correctly, arround 50% for the second and arround 40% for the third. 

![MultiQC](images/bowtie2_se_plot.png)

### STAR

[STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) software is used to aligned reads to a reference genome. 

![MultiQC - Star stats plot](images/star_alignment_plot.png)

- **Uniquely mapped** : successfully aligned reads to one single locus.  
- **Mapped to multiple loci** : accepted reads aligned on 2 to 10 loci.
- **Mapped too many loci** : reads aligned to more than 10 loci. 
- **Unmapped: too many mismatches** : reads having an alignement of too poor quality
- **Unmapped: too short** : less than 66% of reads length (R1+R2) are correctly aligned on the genome. 
- **Unmapped other: other** : other reasons than "too short" or "too many" like for example due to a foreign genome contamination or if reads came from a
from a higly repeated region. 

Alignement is made on correctly barcoded reads.  
Reads passing alignment are uniquely mapped or mapped to multiple loci.   
If this two made up to 60%, the alignemnt is a success. If they made betwwen 60 and 40%, be carrefull and if they made less than 40% the sample may have a problem.  

### Alignment Scores'

Summary of alignment scores accross the pipeline. 

![MultiQC](images/scChIPseq_alignments_plot.png)

- **Deduplicated** : Correct unique reads.
- **Window duplicates** : Remove duplicates by window (if R2 is unmapped).
- **RT duplicates** : Removed RT duplicates which are reads having similar barcode, R2 position but a different R1 position.
- **PCR duplicates** : Removed PCR duplicates which correspond to reads having exactly the same barcode, R1 position and R2 position on the genome.
- **Uniquely mapped not barcoded** : Correctly aligned reads but having an undetermined barcode.
- **Mapped to multiple loci** : Reads aligned on 2 to 10 loci.
- **Unmapped** : All unmaped reads (too many mismatches + too short + other on STAR results).

1) Deduplicated reads   
<span style="color: green">More than 10% : Success</span>  
<span style="color: yellow"> 5-10% : Warning </span>  
<span style="color: red">Less than 5%: Danger</span>  
  
2) RT duplicates  
<span style="color: green">Less than 40% : Success</span>  
<span style="color: yellow"> 40-50% : Warning</span>  
<span style="color: red">More than 50% : Danger</span>  
  
3) PCR duplicates  
<span style="color: green">Less than 40% : Success</span>  
<span style="color: yellow">40-50% : Warning</span>  
<span style="color: red">More than 50% : Danger</span>  
  
4) Window duplicates  
<span style="color: green;">Less than 10% : Success</span>  
<span style="color: yellow;"> 10-15% : Warning</span>  
<span style="color: red;">More than 15% : Danger</span>  


## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `multiqc`**

* `multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser.
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline.

For more information about how to use MultiQC reports, see http://multiqc.info.
