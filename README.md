
# SRDRE (Spatially Resolved Detection of RNA Editing)

The framework SRDRE (Spatially Resolved Detection of RNA Editing) was developed to identify RNA editing events in spatial RNA sequencing data. To eliminate potential genomic single nucleotide polymorphisms (SNPs), a curated catalog of known A-to-I RNA editing events was utilized. For each spot, RNA editing detection was performed using the mapped BAM (Binary Alignment Map) file, which includes aligned sequencing reads and spatial coordinates. The analysis excluded PCR duplicate reads and reads aligned to multiple loci. In summary, SRDRE measured the frequency of A-to-G mismatches at each RNA editing site across spots, enabling in situ visualization of RNA editing events.


## Requirements

- Samtools

- Perl

### Data Preparation

SRD-RE requires the following input data types:

- **Aligned sequencing reads** in sorted BAM format (`input.bam`).

- **Known RNA editing sites** in tab-delimited format , including at least three columns: *Chromosome*, *Coordinate*, *Reference Base* (`knownSite`).

  | Chromosome | Coordinate | Reference Base |
  | ---------- | ---------- | -------------- |
  | chr1       | 1234567    | A              |
  | chr2       | 7654321    | C              |
  | chrX       | 2345678    | G              |


