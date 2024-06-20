##This is the code used to detect somatic mutations between two modules of *Orbicella faveolata*


'''
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=400gb
#SBATCH --job-name=Ofav88
#SBATCH --account=open

source ~/.bashrc

conda activate gatk 

cd ~/scratch/ofav_comparisons

parallel -j 30 gatk Mutect2 -R ~/scratch/ofav_comparisons/um_ofav_v1_softmask.fa \
 -I /storage/group/dut374/default/trinity/bam/20033_SoGv.bam \
 -I /storage/group/dut374/default/trinity/bam/20021_SoGv.bam \
 -tumor 20021 \
 -normal 20033 \
 -O ofav88_{.}.vcf.gz \
 -L {} ::: *.intervals 

 '''
