Bin=`pwd`/../../StereoSeq
samtools=`pwd`/../../tools/samtools

indir=`pwd`/Data
outdir1=`pwd`/Output/1.rmDup
outdir2=`pwd`/Output/2.REs
outdir3=`pwd`/Output/3.REsInfo

echo -e "rmDup\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
mkdir -p $outdir1 && \
perl $Bin/rmDupStStereoSeq.pl --inBam $indir/macaque_cortex_14r.bam --outBam $outdir1/macaque_cortex_14r.bam --samtools $samtools && \
echo -e "rmDup\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

echo -e "REcalling\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
perl $Bin/REcallingStStereoSeq.pl --dataset $indir/Sites_Macaque_macFas5.gz --bam $outdir1/macaque_cortex_14r.bam --outdir $outdir2 --samtools $samtools && \
echo -e "REcalling\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

echo -e "spatialAssign\tStart: $(date +"%H:%M:%S on %d %b %Y")" && \
perl $Bin/spatialAssign.pl --annotation $indir/macaque_cortex_14r.annotation.tsv.gz --input $outdir2/macaque_cortex_14r.REs.gz --outdir $outdir3 --suffix REs.gz && \
echo -e "spatialAssign\tEnd: $(date +"%H:%M:%S on %d %b %Y")"

