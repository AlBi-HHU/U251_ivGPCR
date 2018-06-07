#!/bin/bash
# properties = {properties}

#loads python 3.6.5
module load Python/3.6.5
module load Miniconda3
module load heinz

function handle_qdel {{
    echo "Killed by scheduler"
    mkdir -p $(dirname {jobfailed})
    touch {jobfailed}
}}

trap handle_qdel USR2
trap handle_qdel TERM

{exec_job}

# if the job succeeds, snakemake
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the$
# an error code of 100 is needed since UGER only prevents execution of d$
# job exits with error code 100
cat {jobfinished} &>/dev/null && exit 0 || exit 1

