## Demo Visium

Users can directly repeat all the analysis by following our step-by-step workflow in this file.<br>

We reran this workflow, which took approximately 160 seconds using 1 CPU and 2 GB of memory.



## Data preparation

All files required for detecting RNA A-to-I editing events in Visium data are provided in the `Demo/Visium` directory. Users only need to copy this directory as a whole.

<mark>The only modification to note is that the BAM file should be download directly via link [mouse_brain.bam](https://github.com/JinRcn/SRDRE/raw/refs/heads/main/Demo/Visium/Data/mouse_brain.bam?download=) or by `git lfs pull` .</mark>

Afterward,  replace the LFS pointer file in `Demo/Visium/Data`  and verify the file integrity using `md5sum`.

![DATA](https://github.com/user-attachments/assets/048f14be-d466-4fff-b819-18f8817edc9d)

## Start analysis

Users can directly run the demo with:

```shell
cd Demo/Visium/
# start
sh shell.sh
```

Below we provide a step-by-step explanation of what the script does:

- **Step 0**: Define tool paths and input/output directories

  ```shell
  Bin=`pwd`/../../Visium
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
  perl $Bin/rmDupStVisiumIllumina.pl --inBam $indir/mouse_brain.bam --outBam $outdir1/mouse_brain.bam --samtools $samtools && \
  echo -e "rmDup\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

- **Step 2**: Detection of RNA editing (time-consuming).

  ```shell
  echo -e "REcalling\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
  perl $Bin/REcallingStVisiumIllumina.pl --barcode2slide $indir/visium-v1_coordinates.txt --dataset $indir/REDIportalV2.0_Mouse_mm10.txt.gz --bam $outdir1/mouse_brain.bam --outdir $outdir2 --samtools $samtools && \
  echo -e "REcalling\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

- **Step 3**: Assign spatial annotations to RNA editing sites and calculate RNA editing index

  ```shell
  echo -e "spatialAssign\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
  perl $Bin/spatialAssign.pl --barcode2slide $indir/visium-v1_coordinates.txt --annotation $indir/mouse_brain.annotation.tsv --input $outdir2/mouse_brain.REs.gz --outdir $outdir3 --suffix REs.gz && \
  echo -e "spatialAssign\tEnd: $(date +"%H:%M:%S on %d %b %Y")"
  ```

  

## Final result

![OUTPUT](https://github.com/user-attachments/assets/bc03a5ad-9d62-49cf-96f9-314835af55d2)



All output files can be found in the `Output` directory. The final processed results are in `Output/3.REsInfo`.

- A matrix contains assigned spatial annotations (`mouse_brain.REs.gz`)

  | Coordinate | Chr   | Pos       | RefBase | Coverage | Edited | Area        |
  | :--------- | :---- | :-------- | :------ | :------- | :----- | :---------- |
  | 100,20     | chr1  | 143749562 | A       | 1        | 0      | Isocortex_1 |
  | 100,20     | chr11 | 72209575  | T       | 1        | 0      | Isocortex_1 |
  | 100,20     | chr12 | 100207186 | A       | 5        | 0      | Isocortex_1 |

- A summary table presenting the RNA editing index for each brain area and cell type (`mouse_brain.REI.tsv`).

  | Region       | Coverage | Edited | REI     |
  | :----------- | :------- | :----- | :------ |
  | CA1_CA2      | 16120    | 771    | 0.04783 |
  | CA3          | 14817    | 914    | 0.06169 |
  | Hypothalamus | 60348    | 3734   | 0.06187 |



