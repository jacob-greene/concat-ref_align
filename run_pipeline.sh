#!/bin/sh
#before running snakemake, do in tmux terminal:
module load snakemake/5.19.2-foss-2019b-Python-3.7.4
module load BWA/0.7.17-GCCcore-11.2.0
module load SAMtools/1.16.1-GCC-11.2.0
module load picard/2.25.1-Java-11
module load Java/11.0.2
module load GATK/4.1.4.1-GCCcore-8.3.0-Java-11

#command to run snakemake (remove -np at end when done validating):
snakemake -s subtract_mouse_and_realign.snakefile --latency-wait 60 \
 --keep-going --cluster-config config/cluster_slurm.yamodule load \
  --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
