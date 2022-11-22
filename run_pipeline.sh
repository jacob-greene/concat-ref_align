#!/bin/sh
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml BWA/0.7.17-foss-2018b
ml SAMtools/1.10-GCCcore-8.3.0
ml picard/2.18.29-Java
ml Java/11.0.2
ml GATK/4.1.4.1-GCCcore-8.3.0-Java-11
ml R/3.6.2-foss-2019b-fh1


#command to run snakemake (remove -np at end when done validating):
snakemake -s subtract_mouse_and_realign.snakefile --latency-wait 60 \
 --keep-going --cluster-config config/cluster_slurm.yaml \
  --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
