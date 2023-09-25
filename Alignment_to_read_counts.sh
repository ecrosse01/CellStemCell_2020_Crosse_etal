## Code for the alignment of spatial transcriptome data and read counting

mkdir GENOME
cd GENOME
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gunzip Homo_sapiens.GRCh38.86.gtf.gz

####[Indexing.sh] (Indexing the STAR genome)
#!/bin/sh

#$ -cwd
#$ -pe sharedmem 12
#$ -l h_vmem=8G


. /etc/profile.d/modules.sh

module load igmm/apps/STAR/2.5.1b

ulimit -v

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir dir/GENOME --genomeFastaFiles dir/GENOME/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile dir/GENOME/Homo_sapiens.GRCh38.86.gtf --sjdbOverhang 74

####Trimming###
mkdir output_trimmed

#!/bin/sh
. /etc/profile.d/modules.sh
module load igmm/apps/flexbar/2.5

cd 
for d in dir
do
( cd "$d" && flexbar -r *_1.fastq -p *_2.fastq -a dir/RNA_Refs_Dir/illumina_multiplex.fa -ao 7 --min-read-length 36 --target dir/output_trimmed/"$d" -n 4 )

done

#######Alignment############

#!/bin/sh
#$ -cwd
#$ -pe sharedmem 10
#$ -l h_vmem=8G

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/applications/modulefiles/Community

module load igmm/apps/STAR/2.5.1b

ulimit -v

cd /directorypath/output_trimmed

for read1 in *_1.fastq

do read2=$(echo $read1| sed 's/_1.fastq/_2.fastq/') && STAR --runThreadN 10 --genomeDir dir/GENOME/ --readFilesIn $read1 $read2 --outFileNamePrefix dir/Aligned_reads/$read1
done

##########
###SAM to BAM##


#!/bin/sh

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/applications/modulefiles/Community

module load igmm/apps/samtools/1.6

ulimit -v

cd Aligned_reads

for i in *.sam

do samtools view -b -S -o $i.bam $i

done

####SortBAM###

#!/bin/sh

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/applications/modulefiles/Community

module load igmm/apps/samtools/1.6

ulimit -v

cd /exports/eddie/scratch/s1532897/Aligned_reads

for i in *.bam

do samtools sort $i -o $i.sorted

done


####Indexing###

#!/bin/sh

#$ -cwd
#$ -pe sharedmem 10
#$ -l h_vmem=8G
. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/applications/modulefiles/Community

module load igmm/apps/samtools/1.6

cd dir/Aligned_reads

for i in *.sorted

do samtools index $i

done

###Count the features###  

### BEDCount.sh
wget ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz
gunzip Homo_sapiens.GRCh38.95.gff3.gz
awk '$3 == "gene"' Homo_sapiens.GRCh38.95.gff3 > Homo_sapiens.GRCh38.95.genes.gff3


#!/bin/sh

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/applications/modulefiles/Community

module load igmm/apps/BEDTools/2.27.1

cd /dir/Aligned_reads

for i in *.sorted

do

bedtools multicov -D -p -bams $i -bed dir/Anno/Homo_sapiens.GRCh38.95.gff3> $i.out.coverage

done

####Select the appropriate information i.e. transcripts and counts

#!/bin/sh

for i in *.out.coverage

do

cut -f 9,10 $i | awk -F";" '$1=$1' OFS="\\t" | cut -f 1,8 | sed 's/^.*ID=gene://' > $i.txt

done


cut -f 9,10 CS16V5_1.fastqAligned.out.sam.bam.sorted.out.coverage | awk -F";" '$1=$1' OFS="\\t" | cut -f 1,2,9 | sed 's/^.*Parent=transcript://' | awk -F"=" '$1=$1' OFS="\\t" | cut -f 1,3,4 | awk '$1 ~/ENST/' > CS16V5.txt
