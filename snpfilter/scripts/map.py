map_shell_string = """
TAG=$1
REF=$2
R1=$3
R2=$4
CPU=$5

bwa mem -o $TAG.sam -t $CPU $REF $R1 $R2
samtools faidx $REF
samtools view -bS -@ $CPU -o $TAG.bam $TAG.sam
samtools fixmate -@ $CPU -r -m $TAG.bam $TAG.fm.bam
samtools sort -@ $CPU -o $TAG.st.fm.bam $TAG.fm.bam
samtools markdup -@ $CPU -r $TAG.st.fm.bam $TAG.ud.st.fm.bam
# samtools mpileup -I -t DP,AD -q 20 -Q 20 -g -o $TAG.bcf -f $REF $TAG.ud.st.fm.bam
samtools mpileup -I -t DP,AD -q %d -Q %d -g -o $TAG.bcf -f $REF $TAG.ud.st.fm.bam
samtools index $TAG.ud.st.fm.bam
bcftools index $TAG.bcf

rm $TAG.sam $TAG.bam $TAG.fm.bam $TAG.st.fm.bam
"""