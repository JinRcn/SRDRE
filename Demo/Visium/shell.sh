Bin=`pwd`/../../Visium
samtools=`pwd`/../../tools/samtools

indir=`pwd`/Data
outdir1=`pwd`/Output/1.rmDup
outdir2=`pwd`/Output/2.REs
outdir3=`pwd`/Output/3.REsInfo

echo -e "rmDup\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
mkdir -p $outdir1 && \
perl $Bin/rmDupStVisiumIllumina.pl --inBam $indir/mouse_brain.bam --outBam $outdir1/mouse_brain.bam --samtools $samtools && \
echo -e "rmDup\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

echo -e "REcalling\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
perl $Bin/REcallingStVisiumIllumina.pl --barcode2slide $indir/visium-v1_coordinates.txt --dataset $indir/REDIportalV2.0_Mouse_mm10.txt.gz --bam $outdir1/mouse_brain.bam --outdir $outdir2 --samtools $samtools && \
echo -e "REcalling\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

echo -e "spatialAssign\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
perl $Bin/spatialAssign.pl --barcode2slide $indir/visium-v1_coordinates.txt --annotation $indir/mouse_brain.annotation.tsv --input $outdir2/mouse_brain.REs.gz --outdir $outdir3 --suffix REs.gz && \
echo -e "spatialAssign\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

