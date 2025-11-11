# Demo Stereo-seq

Users can directly repeat all the analysis by following our step-by-step workflow in this file.<br>
We reran this workflow, which took approximately 104 seconds using 1 CPU and 2 GB of memory.



## Data preparation

All files required for detecting RNA A-to-I editing events in Stereo-seq data are provided in the `Demo/StereoSeq` directory. Users only need to copy this directory as a whole.

<mark>The input BAM file are managed by Git LFS and must be downloaded separately, as it will not be downloaded when cloning the repository or downloading it as a ZIP.</mark>

So **the only modification** to note is that the BAM file should be download directly via link [macaque_cortex_14r.bam](https://github.com/JinRcn/SRDRE/raw/refs/heads/main/Demo/StereoSeq/Data/macaque_cortex_14r.bam?download=) or by `git lfs pull` . 

Afterward,  replace the LFS pointer file in `Demo/StereoSeq/Data` and verify the file integrity using `md5sum`.

![DATA](https://github.com/user-attachments/assets/c50c146b-348b-44e7-afcb-def0d509eb39)

## Start analysis

Users can directly run the demo with:

```shell
cd Demo/StereoSeq/
# start
sh shell.sh
```

Below we provide a step-by-step explanation of what the script does:

- **Step 0**: Define tool paths and input/output directories

  ```shell
  Bin=`pwd`/../../StereoSeq
  samtools=`pwd`/../../tools/samtools
  
  indir=`pwd`/Data
  outdir1=`pwd`/Output/1.rmDup
  outdir2=`pwd`/Output/2.REs
  outdir3=`pwd`/Output/3.REsInfo
  ```

- **Step 1**: Remove PCR duplicate reads and reads aligned to multiple loci.

  ```SHELL
  echo -e "rmDup\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
  mkdir -p $outdir1 && \
  perl $Bin/rmDupStStereoSeq.pl --inBam $indir/macaque_cortex_14r.bam --outBam $outdir1/macaque_cortex_14r.bam --samtools $samtools && \
  echo -e "rmDup\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

- **Step 2**: Supervised detection of RNA editing (time-consuming).

  ```shell
  echo -e "REcalling\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
  perl $Bin/REcallingStStereoSeq.pl --dataset $indir/Sites_Macaque_macFas5.gz --bam $outdir1/macaque_cortex_14r.bam --outdir $outdir2 --samtools $samtools && \
  echo -e "REcalling\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

- **Step 3**: Assign spatial annotations to RNA editing sites and calculate RNA editing index

  ```shell
  echo -e "spatialAssign\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
  perl $Bin/spatialAssign.pl --annotation $indir/macaque_cortex_14r.annotation.tsv.gz --input $outdir2/macaque_cortex_14r.REs.gz --outdir $outdir3 --suffix REs.gz && \
  echo -e "spatialAssign\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

  

## Final result
![OUTPUT](https://github.com/user-attachments/assets/604fadd5-82b4-4f08-9f7c-298f684284b7)

All output files can be found in the `Output` directory. The final processed results are in `Output/3.REsInfo`.
- A matrix contains assigned spatial annotations (`macaque_cortex_.REs.gz`)

  | Cod         | Chr   | Pos       | RefBase | Coverage | Edited | Area   | Celltype           |
  | :----------- | :----- | :--------- | :------- | :-------- | :------ | :------ | :----------------- |
  | 20002,77405 | chr8  | 131520619 | T       | 1        | 0      | 14r-l1 | Excitatory neurons |
  | 20003,77362 | chr19 | 55074976  | A       | 1        | 0      | 14r-l1 | Excitatory neurons |
  | 20003,77370 | chr20 | 7610364   | A       | 2        | 0      | 14r-l1 | Excitatory neurons |

- A summary table presenting the RNA editing index for each brain area and cell type (`macaque_cortex_.REI.tsv`).

  | Region | Celltype           | Coverage | Edited | REI     |
  | :------ | :------------------ | :-------- | :------ | :------- |
  | 14r-l2 | Excitatory neurons | 285147   | 27096  | 0.09502 |
  | 14r-l2 | Inhibitory neurons | 63694    | 6402   | 0.10051 |
  | 14r-l3 | Excitatory neurons | 754936   | 75086  | 0.09502 |













