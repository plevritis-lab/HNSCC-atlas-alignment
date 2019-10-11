# AbSeq Tools

## About

This repo houses Python scripts that support BD AbSeq and Single-cell Multiplexing Kit analysis. These scripts are:


`create_ref`: generate fasta and gtf files for BD AbSeq sequences and BD Sample Tags


`abseq_counts`:

* for BD AbSeq analysis
* count AbSeq reads per cell in the BAM alignment file from 10X cellranger, collapse AbSeq reads into molecules, and output data files containing AbSeq molecules per cell count
* optional sample multiplexing analysis: count sample tag reads per cell in BAM alignment file, assign each cell to a sample, and split gene counts (both mRNA and AbSeq) per sample

`demux_gene_counts`:

* for BD Single-cell Multiplexing analysis
* count sample tag reads per cell in the BAM
* assign each cell to a sample and split gene counts per sample

​
​
## Release Notes
------
* V0.1: Initial Release (Jun-15-2018)
* V0.2: Updates (Aug-16-2018)
    1. Added option to generate BDAbSeq reference files based on a list of AbSeq barcode IDs (e.g. AHS0001)
    2. Added new warnings and options in demux_gene_counts for BD Single-cell multiplexing analysis (see multiplexing_tools repo).

## Setup
------
Requires Python 2.7

**Installing abseq_tools with downloaded scripts**

**1. Download and unzip scripts**

**2. Install abseq_tools**

    pip install /path/to/abseq_tools   

NOTE: If permission denied error occurs during pip install, try to install in a virtual environment (see below).

**Install abseq_tools in virtual environment**

**1. Download python (e.g. python 2.7.14 in this example):**  
Open terminal or command prompt and type the following commands:

        cd ~
        mkdir tmp
        cd tmp
        wget https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz
        tar zxvf Python-2.7.14.tgz

**2. Build and install python**

        cd Python-2.7.14
        ./configure --prefix=$HOME/opt/python-2.7.14
        make
        make install

**3. Add path**

        export PATH=$HOME/opt/python-2.7.14/bin:$PATH

**4. Create python virtual environment**  
Make sure to change 'username' to your username

        virtualenv -p /home/username/opt/python-2.7.14/bin/python myvenv

**5. Activate virtual environment (this needs to be done before you run the script every time)**

        source myvenv/bin/activate

**6. Install abseq_tools**  
To install through bitbucket:

        pip install git+https://bitbucket.org/CRSwDev/abseq_tools

To install through donwloaded scripts:

        pip install /path/to/abseq_tools    

# Workflow for BD AbSeq Analysis
-----
1. Generate BD AbSeq reference files (.gtf and .fasta) using `AbSeq_tools`
2. Build refernece using 10X Cell Ranger
3. Process RNA and AbSeq sequencing libraries separately: 
    * RNA library: run 10X cellranger pipeline for read alignment and cell-gene count table generation using desired reference genome
    * AbSeq library: run 10X cellranger pipeline for read alignment using BD AbSeq reference 
4. Use `AbSeq_tools` to count AbSeq reads/molecules and generate cell-gene data table

### How to run `create_ref`

`create_ref --bdabseq --reflist=abseq_list.txt`

#### Input
[Required] Specify AbSeq reference generation option.

`--bdabseq`

[Required] One of the following options are required to specify what type of reference to generate.

Option 1: A list of BD AbSeq barcode IDs in a txt file, one barcode ID per line (e.g. AHS0001). For example:

        AHS0001
        AHS0002
        AHS0003
        AHS0004
        AHS0021

If this barcode ID list is supplied, the script will generate a fasta file and a gtf file containing the list of AbSeqs specified.

`--reflist=abseq_list.txt`

Option 2:  FASTA file containing list of AbSeq barcode sequences and custom names for each. If this option is supplied, the script will generate a gtf file.

`--fasta=new_abseq.fasta`


[Optional] Specify the name of the output fasta and gtf files. Default name is BDAbSeq (i.e. BDAbSeq.fasta, BDAbSeq.gtf)

`--genome=newName`

#### Output

* Fasta and gtf files containing BD AbSeq sequences.


### How to run `abseq_counts` 

​
        
    abseq_counts \
        --bam=possorted_genome_bam.bam \
        --filtered_barcodes=./outs/filtered_gene_bc_matrices/GRCh38 \
        --gtf=BDAbSeq.gtf \
        --rna_counts=./outs/filtered_gene_bc_matrices/GRCh38 \
        --output=path/to/output
            
#### Input

[Required] BAM alignment file of the BD AbSeq library from 10X Cell Ranger 

`--bam=possorted_genome_bam.bam`

[Required] Directory containing `barcodes.tsv`, which is the list of cell barcodes of the filtered cells based on the mRNA data. 

`--filtered_barcodes=./outs/filtered_gene_bc_matrices/GRCh38`

[Required] gtf file of BD AbSeq sequences used for alignment.

`gtf=BDAbSeq.gtf`

[Required] Output directory. The script will automatically generate the specified directory with output files. Must be a new directory.

`--output=/path/to/output`

[Optional] Directory containing RNA gene cell matrix files from 10X Cell Ranger (`matrix.mtx`, `barcodes.tsv`, `genes.tsv`). If supplied, the script will combined the mRNA and newly generated AbSeq data table.

`rna_counts=./outs/filtered_gene_bc_matrices/GRCh38`

[Optional] Molecule info file containig RNA molecule info from 10X Cell Ranger (`molecule_info.h5`). If supplied, the script will generate a CSV file containing the total AbSeq and mRNA read and molecule counts per cell barcode prior to cell label filtering (`Total_RnaAbSeqPerCell.csv`).

`rna_info=./outs/molecule_info.h5`

[Optional] Save raw molecule counts of each AbSeq (columns) per cell barcode (rows) as a CSV (`AbSeq_MolsPerCell.csv`)

`--save_raw_counts`

[Optional] For AbSeq + sample multiplexing analysis only. BAM alignment file from a single-cell experiment where cells are identified by a CB tag OR csv file containing sample tag counts for each cell label. 

`--demux_counts_bam=./sampletag/outs/possorted_genome_bam.bam`


#### Output
* Molecule Gene cell matrix files that contain only AbSeq counts (`AbSeq_MolsPerCell` folder with `matrix.mtx`, `barcodes.tsv`, `genes.tsv`)
* Molecule Gene cell data table that contains only AbSeq counts (`AbSeq_filteredMolsPerCell.csv`)
* Histogram of AbSeq counts 
* Metrics file containing metrics for each AbSeq barcode (`stats.csv`)
* [Optional: output of `--rna_counts`] Molecule Gene cell matrix files that combines RNA and AbSeq counts (`AbSeqRNA_MolsPerCell` folder with `matrix.mtx`, `barcodes.tsv`, `genes.tsv`)
* [Optional: output of `--rna_info`] Molecule Gene cell data table that contains only AbSeq counts from all cell barcodes (`AbSeq_MolsPerCell.csv`)
* [Optional: output of `--save_raw_counts`] `AbSeq_MolsPerCell.csv` containing molecule counts of each AbSeq (columns) per cell barcode (rows) prior to cell barcode filtering
* [Optional: output of `--demux_counts_bam`] If optional sample multiplexing analysis is selected, molecule gene cell matrix files will be split by BD sample tags

​
​
## Workflow for BD Single-cell Multiplexing Analysis
------
1. Generate Sample Tag reference files
2. Combine Sample Tag reference with your desired reference genome
3. Run single-cell RNASeq pipeline for read alignment and cell-gene count table generation
4. Analyze Sample Tag reads and assign cells to samples

### How to run `demux_gene_counts`

​
        `demux_gene_counts \
            --bam=possorted_genome_bam.bam \
            --counts=./outs/filtered_gene_bc_matrices/GRCh38 \
            --output=path/to/output`

##### Input

[Required] Sample tag count information (one of these files are required): 
BAM alignment file from a single-cell experiment where cells are identified by a CB tag OR csv file containing sample tag counts for each cell label. 

`--bam=possorted_genome_bam.bam`    OR    `--tag_counts=SampleTag_ReadsPerCell.csv`


[Required] List of cell barcodes of the filtered cells based on the mRNA data It can be provided in one of the following ways:

* as a directory containing mRNA gene cell matrix files (matrix.mtx, barcodes.tsv, and genes.tsv files, or an h5 file, or a Rhapsody Expression_Data file, or a CSV file with Genes in columns and cells in rows.)
 
    `--counts=path/to/expression_matrix.mtx`    
​    

* as a list of cell barcodes in a CSV or TXT file

    `--barcodes=barcodes.csv`
​    

[Required] Output directory. The script will automatically generate the specified directory with output files. Must be a new directory.

`--output=/path/to/output`

[Optional] GTF file containing sample tag references. This is only needed if a sample tag reference other than the default reference files (e.g.: 'BDSampleTags.fasta', 'BDSampleTags.gtf' 'BDSampleTagsMM.fasta', 'BDSampleTagsMM.gtf') were used to generate BAM alignment results.   

`--tag_gtf=new_sample_tag.gtf`

[Optional] Save raw read counts of each sample tag (columns) per cell barcode (rows) as a CSV

`--save_raw_counts`

[Optional] Specify a threshold for minimum number of reads per Sample Tag per cell. Counts that do not reach this threshold will be discarded.  Used for noise reduction algorithm to demultiplex sample tags.

`--min_cell_sample_read=20`


##### Output 

* `cell2Label.csv` containing read counts per cell for each BD Sample Tag and sample tag assignment
* `Sample_Tag_Calls.csv` containing sample tag assignment of each cell barcode
* `Stats.csv` containing statistics of read counts, cell numbers, and signal to noise ratio metrics in each BD Sample Tag
* [Optional: output of `--counts`] Sub-folders that splits the original Molecule gene cell matrix into multiple Molecule Gene cell matrix file by BD Sample Tag
* [Optional: output of `--save_raw_counts`] tagReadCounts.csv containing raw read counts of each sample tag (columns) per cell barcode (rows)


### How to run `ceate_ref`:

`create_ref --bdgenome --genome=BDTags`

##### Input

[Required] Specify whether to generate reference files for human (--bdgenome) or mus musculus (--bdgenome_mm) sample tags

`--bdgenome` OR `-bdgenome_mm`

[Optional] Specify the name of the output fasta and gtf files. Default name is BDSampleTags (i.e. BDSampleTags.fasta, BDSampleTags.gtf)

`--genome=BDTags`


##### Output:

* Fasta and gtf files containing BD Sample Tag sequences.

## Example of BD AbSeq Workflow
-----

### Using 10x Cell Ranger Gene-Barcode Matrices with AbSeq

1. **Generate reference files for BD AbSeq (BDAbSeq.fasta and BDAbSeq.gtf) and build a referece**

    * Prepare a list of BD AbSeq barcode IDs in a txt file (e.g. reflist.txt), one barcode ID (e.g. AHS0001) per line 
    * Generate reference files using `create_ref` in `AbSeq_tools`
    
            create_ref \
                --bdabseq 
                --reflist=relist.txt
                    ​
  
    * Build BD AbSeq reference using Cell Ranger
      See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references> 

                cellranger mkref \
                    --genome=BDAbSeq \
                    --fasta=BDAbSeq.fasta \
                    --genes=BDAbSeq.gtf  
                    ​
  
  
2. **Run Cell Ranger for RNA and AbSeq library separately**

    **Generate fastq files for RNA and AbSeq library separately:**
    
            cellranger mkfastq \
                --id=${MKFASTQ_OUTDIR} \
                --run=${BCL_FOLDER} \
                --samplesheet=${SAMPLE_SHEET}  
                ​
        
    **Run cellranger count for RNA library:**
            
            cellranger count \
                --id=${COUNT_OUTDIR_RNA} \
                --transcriptome=GRCh38 \
                --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                --sample=myRNAsample  
                ​
        
    **Run cellranger count for AbSeq library specifying trimmed r2-length:**
    
    See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count>
    
            cellranger count \
                --id=${COUNT_OUTDIR_AbSeq} \
                --transcriptome=BDAbSeq \
                --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                --sample=myAbSeqsample
                --r2-length=60  
                ​
                
    **NOTE: If RNA and AbSeq libraries were prepared with the same index:**
            
    **Run all fastq files through Cell Ranger twice, once aligning with desired reference genome, and then repeat with AbSeq reference:**

                cellranger count \
                    --id=${COUNT_OUTDIR_RNA} \
                    --transcriptome=GRCh38 \
                    --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                    --sample=myRNA_AND_AbSeq_sample  
                ​
                cellranger count \
                    --id=${COUNT_OUTDIR_AbSeq} \
                    --transcriptome=BDAbSeq \
                    --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                    --sample=myRNA_AND_AbSeq_sample
                    --r2-length=60  
                ​
    


3. **Run abseq_counts**
    **Once alignment is completed with Cell Ranger, run `abseq_counts` using Cell Ranger output files to generate AbSeq cell-gene data table**

            abseq_counts \
                --bam=${COUNT_OUTDIR_AbSeq}/outs/possorted_genome_bam.bam \
                --rna_counts=${COUNT_OUTDIR_RNA}/outs/filtered_gene_bc_matrices/GRCh38/ \
                --filtered_barcodes=${COUNT_OUTDIR_RNA}/outs/filtered_gene_bc_matrices/GRCh38/ \
                --gtf=path/to/gtf-file/BDAbSeq.gtf \
                --output=path/to/output  
                ​


### Using 10x Cell Ranger Gene-Barcode Matrices with AbSeq and BD Sample Tags

1. **Generate reference files for BD AbSeq (BDAbSeq.fasta and BDAbSeq.gtf) and build a referece**

    * Prepare a list of BD AbSeq barcode IDs in a txt file (e.g. reflist.txt), one barcode ID (e.g. AHS0001) per line 
    * Generate reference files using `create_ref` in `AbSeq_tools`
    
            create_ref \
                --bdabseq 
                --reflist=relist.txt
                    ​
  
    * Build BD AbSeq reference using Cell Ranger
      See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references> 

                cellranger mkref \
                    --genome=BDAbSeq \
                    --fasta=BDAbSeq.fasta \
                    --genes=BDAbSeq.gtf  
                    ​
2. **Generate reference files for BD Sample Tags (BDSampleTags.fasta, BDSampleTags.gtf) using `create_ref` in `AbSeq_tools`**

        create_ref --bdgenome
        ​
  
3. **Download 10x reference**

    See <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>


4. **Generate a reference that includes Sample Tag sequences**

    See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references>

        cellranger mkref \
            --genome=GRCh38  \
            --fasta=${10X_REFERENCE}/fasta/genome.fa \
            --genes=${10X_REFERENCE}/genes/genes.gtf \
            --genome=BDSampleTags \
            --fasta=BDSampleTags.fasta \
            --genes=BDSampleTags.gtf  
            ​

5. **Run Cell Ranger for mRNA+BDSampleTag and AbSeq library separately**

    **Generate fastq for RNA+BDSampleTags and AbSeq library:**
    
        cellranger mkfastq \
            --id=${MKFASTQ_OUTDIR} \
            --run=${BCL_FOLDER} \
            --samplesheet=${SAMPLE_SHEET}  
            ​
    
    **Run cellranger count for RNA+BDSampleTag library:**

        cellranger count \
            --id=${COUNT_OUTDIR_RNA} \
            --transcriptome=GRCh38_and_BDSampleTags \
            --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
            --sample=myRNAandBDSampleTagsample  
            ​
            
            
    **Run cellranger count for AbSeq library specifying trimmed r2-length:**
    
    See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count>
    
        cellranger count \
            --id=${COUNT_OUTDIR_AbSeq} \
            --transcriptome=BDAbSeq \
            --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
            --sample=myAbSeqsample
            --r2-length=60  
            ​


6. **Run abseq_counts with option to demultiplex samples**

        abseq_counts \
            --bam=${COUNT_OUTDIR_AbSeq}/outs/possorted_genome_bam.bam \
            --rna_counts=${COUNT_OUTDIR_RNA}/outs/filtered_gene_bc_matrices/GRCh38/ \
            --filtered_barcodes=${COUNT_OUTDIR_RNA}/outs/filtered_gene_bc_matrices/GRCh38/ \
            --gtf=path/to/gtf-file/BDAbSeq.gtf \
            --output=path/to/output \
            --demux_counts_bam=${COUNT_OUTDIR_RNA}/outs/possorted_genome_bam.bam  
            ​
            


## Examples of BD Single-Cell Multiplexing Kit Analysis
------

### Using 10x Cell Ranger Gene-Barcode Matrices

1. **Generate reference files for BD Sample Tags (BDSampleTags.fasta, BDSampleTags.gtf) using `create_ref` in `AbSeq_tools`**

        create_ref --bdgenome
        ​

2. **Download 10x reference**

    See <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>

3. **Generate a reference that includes Sample Tag sequences**
    
    * If mRNA and sample tag reads are sequenced in the SAME library, generate a combined reference that includes both human genome and sample tag reference

        See <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references>

            cellranger mkref \
                --genome=GRCh38  \
                --fasta=${10X_REFERENCE}/fasta/genome.fa \
                --genes=${10X_REFERENCE}/genes/genes.gtf \
                --genome=BDSampleTags \
                --fasta=BDSampleTags.fasta \
                --genes=BDSampleTags.gtf
       ​
       
    * If mRNA and sample tag reads are sequenced in SEPARATE libraries, generate 2 separate references, one containing human genome reference, and a separate one containing sample tag reference
    ​


            cellranger mkref \
                --genome=GRCh38  \
                --fasta=${10X_REFERENCE}/fasta/genome.fa \
                --genes=${10X_REFERENCE}/genes/genes.gtf \

            cellranger mkref 
                --genome=BDSampleTags \
                --fasta=BDSampleTags.fasta \
                --genes=BDSampleTags.gtf
                ​
                
4. **Run Cell Ranger to generate mRNA data table and sample tag alignment BAM file**

        cellranger mkfastq \
            --id=${MKFASTQ_OUTDIR} \
            --run=${BCL_FOLDER} \
            --samplesheet=${SAMPLE_SHEET}
            ​
                
    * If mRNA and sample tag reads are sequenced in the SAME libraries, run CellRanger on ALL fastq files together:
            
            cellranger count \
                --id=${COUNT_OUTDIR} \
                --transcriptome=GRCh38_and_BDSampleTags \
                --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                --sample=mysample
                ​
                        
    * If mRNA and sample tag reads are sequenced in SEPARATE libraries, run CellRanger on each set of fastq files separately:
             
            cellranger count \
                --id=${COUNT_OUTDIR} \
                --transcriptome=GRCh38 \
                --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                --sample=myRNASample
                
            cellranger count \
                --id=${COUNT_OUTDIR} \
                --transcriptome=BDSampleTags \
                --fastqs=${MKFASTQ_OUTDIR}/outs/fastq_path \
                --sample=myBDSampleTagSample
                        ​
                        
5. **Run demux_gene_counts**

        demux_gene_counts \
            --bam=${COUNT_OUTDIR}/outs/possorted_genome_bam.bam \
            --counts=${COUNT_OUTDIR}/outs/filtered_gene_bc_matrices/GRCh38 \
            --output=path/to/output


### Using Sample Tag Read Counts Per Cell data table

1. **Run demux_gene_counts**

        demux_gene_counts --tag_counts=SampleTag_ReadsPerCell.csv --barcodes=barcodes.csv --output=/path/to/output


## Loading Data into DataView 
------
* For BD Single Cell Multiplexing analysis:
    * Load filtered gene cell matrix output from Cell Ranger (filtered_gene_bc_matrices -> matrix.mtx, genes.tsv, barcodes.tsv) into DataView.
        * For more information, please referece the BD Data View chapter -> Load data section in [BD single cell genomics bioinformatics handbook](http://www.bd.com/documents/guides/user-guides/GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf)
    * Load Sample_Tag_Calls.csv output from demux_gene_counts to annotate cells by BD Sample Tags
        * For more information, please referece the BD Data View chapter -> Analyzing a multiplexed sample with BD Data View in [BD single cell genomics bioinformatics handbook](http://www.bd.com/documents/guides/user-guides/GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf)

* For BD AbSeq Samples:
    * Load gene cell matrix output from `AbSeq_tools` (AbSeqRNA_MolsPerCell -> matrix.mtx, genes.tsv, barcodes.tsv) into DataView. This data matrix contains both AbSeq and RNA molecule counts.
    * To selectively analyze RNA and AbSeq markers, use the automatically generated AbSeq gene-set and 'Filter data table' function.
        * For more information, please referece the BD Data View chapter -> Filter data table in [BD single cell genomics bioinformatics handbook](http://www.bd.com/documents/guides/user-guides/GMX_BD-Rhapsody-genomics-informatics_UG_EN.pdf)

## Troubleshooting
------
#### Common Error Messages
 ------

`Could not find any reads mapped to sample tags. Please check that the correct gtf info is supplied.`
 
 ------
The script looks for default `seqname (name of chromosome or scaffold)` `BDSampleTags` for human sample tags and `BDSampleTagsMM` for mouse sample tags, assuming that the default sample tag reference files were not used (generated by `create_ref`).
Potential causes of the error:

 1. If a custom gtf file with a different seqname is used, please make sure to supply a correct gtf file using the input option:`--tag_gtf=custom_gtf_file.gtf`

 2. Check to make sure that the correct BAM file containing sample tag reads is supplied in the input. 
​
       

 ------

`Sample tag read counts are too low for analysis. Please refer to the log file demux_counts.log for details.`
 
 ------
This error message indicates that the sum of each sample tag reads in all cells is less than 1000 for all sample tags. This can occur when the wrong bam file is supplied, or insufficient sequencing of the sample tag library.
If Cell Ranger was used to generate thebam file, please also refer back to the web_summary.html file and make sure that sufficient number of reads (>600 sample tag reads per cell) were mapped to the sample tag library.
​

#### Common Warning Messages
 ------

`Warning: High percentage ({X}%) of cells have no sample tag assignment`
 
 ------
This warning message appears when more than 10% of cells have no sample tag assignment (NoSample). 
The script will still run to completion, but the warning indicates that the noise in sample tags are too high to confidently assign cells to specific sample tags. 
Please check cell2Label.csv and stats.csv for more details.
​

 ------

`Warning: Low signal to noise ratio in Sample Tag {Y} = {X}`

 ------
This warning message appears when the mean signal to noise ratio of sample tag {Y} is less than 10. Low signal to noise ratio can lead to inaccurate sample tag calls. 
The script will still run to completion, but the warning indicates that the noise in sample tags can be too high to confidently assign cells to specific sample tags. 
Please check cell2Label.csv and stats.csv for more details.
​

       



## License

[![Creative Commons License](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).

## Contact

Contact <techsupport@bdgenomics.com> if you have questions

