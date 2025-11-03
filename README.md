
# SRDRE (Spatially Resolved Detection of RNA Editing)

The framework SRDRE (Spatially Resolved Detection of RNA Editing) was developed to identify RNA editing events in spatial RNA sequencing data. To eliminate potential genomic single nucleotide polymorphisms (SNPs), a curated catalog of known A-to-I RNA editing events was utilized. For each spot, RNA editing detection was performed using the mapped BAM (Binary Alignment Map) file, which includes aligned sequencing reads and spatial coordinates. The analysis excluded PCR duplicate reads and reads aligned to multiple loci. In summary, SRDRE measured the frequency of A-to-G mismatches at each RNA editing site across spots, enabling in situ visualization of RNA editing events.

![image](https://github.com/user-attachments/assets/336e7f97-03c7-4623-9b4b-b156d92b2650)

- [Data Preparation](#data-preparation)
- [Stereo-seq](#stereo-seq-data)
- [Visium](#visium-data)

The applicability of our tool has been proven on two widely used spatial transcriptomics platforms, **Stereo-seq** and **Visium.** We uploaded the test data in the [Demo](https://github.com/JinRcn/SRDRE/tree/main/Demo) directory, and users can directly repeat the analysis by following our step-by-step workflow.<br>
**Note**: The BAM files are managed by Git LFS, so please download them directly via link below:  <br>
          [StereoMacaque](https://github.com/JinRcn/SRDRE/raw/refs/heads/main/Demo/StereoSeq/Data/macaque_cortex_14r.bam?download=) and 
          [VisiumMouse](https://github.com/JinRcn/SRDRE/raw/refs/heads/main/Demo/Visium/Data/mouse_brain.bam?download=).

<br>
Use the tool instantly — no installation required.
<br>

## Requirements

- [Samtools](https://github.com/samtools/samtools)

- [Perl](https://www.perl.org)
<br>

### Data Preparation

SRDRE requires the following input data types:

- **Aligned sequencing reads** in sorted BAM format (`input.bam`).

- **Known RNA editing sites** in tab-delimited format , including at least three columns: *Chromosome*, *Coordinate*, *Reference Base* (`knownSite`).

  | Chromosome | Coordinate | Reference Base |
  | :---------- | :---------- | :-------------- |
  | chr1       | 1234567    | A              |
  | chr2       | 7654321    | C              |
  | chrX       | 2345678    | G              |
<br>

### Stereo-seq data

The detection of RNA editing events in Stereo-seq data involves three main steps:

- **Step 1**: Remove PCR duplicate reads and reads aligned to multiple loci.

  **Input** 
  
  - A sorted BAM file aligned to the reference genome (`input.bam`).

  **Command**
  
  ```perl
  perl StereoSeq/rmDupStStereoSeq.pl --inBam <input.bam> --outBam <output.bam> --samtools <samtools>
  ```
  
  | Argument       | Description                                                  |
  | :------------- | :------------------------------------------------------------ |
  | `--inBam `     | The input BAM file. A sorted BAM file with alignments ordered by leftmost coordinates is required as input. |
  | `--outBam `    | The output BAM file after removing duplicates and multi-mapped reads. |
  | `--samtools `  | The path to the `samtools` executable.                       |
  
  **Output**  
  
  - A BAM file after removing PCR duplicates and multi-mapped reads (`sample.bam`).


- **Step 2**: Detection of RNA editing.

  **Input** 
  
  - The BAM file generated from Step 1.
  
  - A known A-to-I RNA editing events dataset.
  
  **Command**
  
  ```perl
  perl StereoSeq/REcallingStStereoSeq.pl --dataset <knownSite> --bam <sample.bam> --outdir <results> --samtools <samtools> --suffix bam --phred 33 --qual_cutoff 20
  ```
  
  | Argument           | Description                                                  |
  | :---------------   | :------------------------------------------------------------ |
  | `--dataset `       | The path to the known A-to-I RNA editing dataset (e.g. Dataset/REDIportalV2.0_Mouse_mm10.txt.gz) |
  | `--bam `           | The input BAM file, generated from Step 1.                   |
  | `--outdir `        | The path to the output directory for results.                |
  | `--samtools `      | The path to the `samtools` executable.                       |
  | `--suffix `        | The suffix of the input file. [Default: bam]                 |
  | `--phred `         | The Phred quality score encoding. [Default: 33]              |
  | `--qual_cutoff `   | The quality cutoff for base calling. [Default: 20]           |
  
  **Output**  
  
  - A count matrix of all detected RNA editing events in the spatial data (`.REs.gz`).
    Each entry contains the X and Y spatial coordinates (**Coordinate**), where **Coverage** indicates the total number of reads covering the locus, while **Edited** represents the number of edited reads at that site. **Chr**, **Pos**, and **RefBase** denote the chromosome identifier, exact genomic position, and reference base of the known A-to-I editing site, respectively.
  
    | Coordinate   | Chr   | Pos       | RefBase | Coverage | Edited |
    | :----------- | :----- | :--------- | :------- | :-------- | :------ |
    | 20002,77405 | chr8  | 131520619 | T       | 1        | 0      |
    | 20003,77362 | chr19 | 55074976  | A       | 1        | 0      |
    | 20003,77370 | chr20 | 7610364   | A       | 2        | 0      |
  
  - A matrix recording the detected base and its Phred quality score for editing site within the given coordinates (`.sam2base.gz`). **CoverBase** records the base detected at the editing site, and **Phred** represents the Phred quality score used for quality control.
  
    | Coordinate    | Chr   | Pos       | RefBase | CoverBase | Phred |
    | :----------- | :----- | :--------- | :------- | :--------- | :----- |
    | 20002,77405 | chr8  | 131520619 | T       | T         | G     |
    | 20003,77362 | chr19 | 55074976  | A       | A         | G     |
    | 20003,77370 | chr20 | 7610364   | A       | AA        | D?    |

​	
- **Step 3**:  Assign spatial annotations to RNA editing sites and calculate RNA editing index

  **Input** 

  - The `.REs.gz` file generated from Step 2.
  
  - An annotation file of the Stereo-seq data.
  
  **Command**
  
  ```perl
  perl StereoSeq/spatialAssign.pl --annotation <annotations.tsv> --input <site.REs.gz> --outdir <results> --suffix REs.gz
  ```
  
  | Argument        | Description                                                  |
  | :-------------- | :------------------------------------------------------------ |
  | `--annotation `  | The annotation file of the Stereo-seq data, including spatial coordinates, cell-type annotations, etc. |
  | `--input `       | The input file containing the detected RNA editing sites, generated from Step 2. |
  | `--outdir `      | The path to the output directory for results.                |
  | `--suffix `      | The suffix of the input file. [Default: REs.gz]              |
  
  **Output**  
  
  - A matrix contains assigned spatial annotations (`.REs.gz`)
  
    | Coordinate    | Chr   | Pos       | RefBase | Coverage | Edited | Area   | Celltype           |
    | :----------- | :----- | :--------- | :------- | :-------- | :------ | :------ | :------------------ |
    | 20002,77405 | chr8  | 131520619 | T       | 1        | 0      | 14r-l1 | Excitatory neurons |
    | 20003,77362 | chr19 | 55074976  | A       | 1        | 0      | 14r-l1 | Excitatory neurons |
    | 20003,77370 | chr20 | 7610364   | A       | 2        | 0      | 14r-l1 | Excitatory neurons |
  
  - A summary table presenting the RNA editing index for each brain area and cell type (`.REI.tsv`).
  
    | Region | Celltype           | Coverage | Edited | REI     |
    | :------ | :------------------ | :-------- | :------ | :------- |
    | 14r-l2 | Excitatory neurons | 285147   | 27096  | 0.09502 |
    | 14r-l2 | Inhibitory neurons | 63694    | 6402   | 0.10051 |
    | 14r-l3 | Excitatory neurons | 754936   | 75086  | 0.09502 |

<br>

### Visium data

The only difference to note in Visium is that a **coordinate file** in tab-delimited format is **required** to map barcodes to their corresponding slide coordinates (`barcodes.tsv`).

| barcode          | x    | y    |
| :---------------- | :---- | :---- |
| AAACAACGAATAGTTC | 17   | 1    |
| AAACAAGTATCTCCCA | 103  | 51   |
| AAACAATCTACTAGCA | 44   | 4    |
| AAACACCAATAACTGC | 20   | 60   |

<br>
The detection of RNA editing events in Visium data follows a process similar to that for Stereo-seq.

- **Step 1**: Remove PCR duplicate reads and reads aligned to multiple loci.

  **Input** 

  - A sorted BAM file aligned to the reference genome (`input.bam`).

  **Command**

  ```perl
  perl Visium/rmDupStVisiumIllumina.pl --inBam <input.bam> --outBam <output.bam> --samtools <samtools>
  ```
  | Argument      | Description                                                  |
  | :------------- | :------------------------------------------------------------ |
  | `--inBam `    | The input BAM file. A sorted BAM file with alignments ordered by leftmost coordinates is required as input. |
  | `--outBam `   | The output BAM file after removing duplicates and multi-mapped reads. |
  | `--samtools ` | The path to the `samtools` executable.                       |
  
  **Output**  

  - A BAM file after removing PCR duplicates and multi-mapped reads (`sample.bam`).


- **Step 2**:  Detection of RNA editing 

  **Input** 

  - The BAM file generated from Step 1.
  
  - A known A-to-I RNA editing events dataset.
  
  - A barcode-to-slide coordinates file.
  
  **Command**
  
  ```perl
  perl Visium/REcallingStVisiumIllumina.pl --barcode2slide <barcodes.tsv> --dataset <knownSite> --bam <sample.bam> --outdir <results> --samtools <samtools> --suffix bam --phred 33 --qual_cutoff 20
  ```
  
  | Argument           | Description                                                  |
  | :----------------- | :------------------------------------------------------------ |
  | `--barcode2slide ` | The file records the correspondence between barcodes and spatial coordinates. |
  | `--dataset `       | The path to the known A-to-I RNA editing dataset (e.g. Dataset/REDIportalV2.0_Mouse_mm10.txt.gz) |
  | `--bam `           | The input BAM file, generated from Step 1.                   |
  | `--outdir `        | The path to the output directory for results.                |
  | `--samtools `      | The path to the `samtools` executable.                       |
  | `--suffix `        | The suffix of the input file. [Default: bam]                 |
  | `--phred `         | The Phred quality score encoding. [Default: 33]              |
  | `--qual_cutoff `   | The quality cutoff for base calling. [Default: 20]           |
  
  **Output**  
  
  - A count matrix of all detected RNA editing events in the spatial data (`mouse_brain.REs.gz`).
    Each entry contains the X and Y spatial coordinates (**Coordinate**), where **Coverage** indicates the total number of reads covering the locus, while **Edited** represents the number of edited reads at that site. **Chr**, **Pos**, and **RefBase** denote the chromosome identifier, exact genomic position, and reference base of the known A-to-I editing site, respectively.

    | Coordinate    | Chr   | Pos       | RefBase | Coverage | Edited |
    | :------ | :----- | :--------- | :------- | :-------- | :------ |
    | 100,20 | chr1  | 143749562 | A       | 1        | 0      |
    | 100,20 | chr11 | 72209575  | T       | 1        | 0      |
    | 100,20 | chr12 | 100207186 | A       | 5        | 0      |
  
  - A matrix recording the detected base and its Phred quality score for editing site within the given coordinates (`mouse_brain.sam2base.gz`).
    **CoverBase** records the base detected at the editing site, and **Phred** represents the Phred quality score used for quality control.  

    | Coordinate    | Chr   | Pos       | RefBase | CoverBase | Phred |
    | :------ | :----- | :--------- | :------- | :--------- | :----- |
    | 100,20 | chr1  | 143749562 | A       | A         | F     |
    | 100,20 | chr11 | 72209575  | T       | T         | F     |
    | 100,20 | chr12 | 100207186 | A       | AAAAA     | FF:FF |


- **Step 3**:  Assign spatial annotations to RNA editing sites and calculate RNA editing index

  **Input** 

  - The `.REs.gz` file generated from Step 2.
  
  - An annotation file of the Visium data.
  
  - A barcode-to-slide coordinates file.
  
  **Command**
  
  ```perl
  perl Visium/spatialAssign.pl --barcode2slide <barcodes.tsv> --annotation <annotations.tsv> --input <site.REs.gz> --outdir <results> --suffix REs.gz
  ```
  | Argument        | Description                                                  |
  | :-------------- | :------------------------------------------------------------ |
  | `--annotation `  | The annotation file of the Stereo-seq data, including spatial coordinates, cell-type annotations, etc. |
  | `--input `       | The input file containing the detected RNA editing sites, generated from Step 2. |
  | `--outdir `      | The path to the output directory for results.                |
  | `--suffix `      | The suffix of the input file. [Default: REs.gz]              |
  
  **Output** 
  
  - A matrix contains assigned spatial annotations (`mouse_brain.REs.gz`).
  
    | Coordinate  | Chr   | Pos       | RefBase | Coverage | Edited | Area        |
    | :------ | :----- | :--------- | :------- | :-------- | :------ | :----------- |
    | 100,20 | chr1  | 143749562 | A       | 1        | 0      | Isocortex_1 |
    | 100,20 | chr11 | 72209575  | T       | 1        | 0      | Isocortex_1 |
    | 100,20 | chr12 | 100207186 | A       | 5        | 0      | Isocortex_1 |
  
  - A summary table presenting the RNA editing index for each brain area and cell type (`mouse_brain.REI.tsv`).
  
    | Region       | Coverage | Edited | REI     |
    | :------------ | :-------- | :------ | :------- |
    | CA1_CA2      | 16120    | 771    | 0.04783 |
    | CA3          | 14817    | 914    | 0.06169 |
    | Hypothalamus | 60348    | 3734   | 0.06187 |

<br>

### **[NOTES]**

Ensure that the paths to the Perl scripts and the `samtools` executable are correct. Replace `input.bam`, `output.bam`, the `knownSite` dataset, `results` path, and other placeholders with actual values specific to your data and environment. 

The `-phred 33` option assumes that the quality scores are encoded using the Illumina 1.8+ format. If your data uses a different encoding, adjust this parameter accordingly. The `-qual_cutoff 20` option sets a quality threshold for base calling, which can be adjusted based on the quality of your sequencing data.

<br>

## Citation

<br>

## Tools comparison
[Comparative analysis](https://github.com/JinRcn/SRDRE/blob/main/Other/comparison.md)
<br>

## Contact

If you have any comments or suggestions, please contact Dr Huang (huangjinrong@genomics.cn).



## License

**SRD-RE is free for academic use only.**
