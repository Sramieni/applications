#!/bin/bash
mv job_* old
for l in `cat lists.txt`;
    do sbatch ~/slurm/saveFiguresToPS.sh $l;
done;
