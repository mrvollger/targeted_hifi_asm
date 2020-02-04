DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -s env.cfg ]; then
	source env.cfg
else 
	source $DIR/env.cfg
fi



SNAKEFILE=$DIR/targeted_hifi_asm.smk.py
JOBNUM=100
WAITTIME=100
LOGDIR=log
RETRY=0
mkdir -p $LOGDIR

snakemake -s $SNAKEFILE \
        --drmaa " -P eichlerlab \
                -q eichler-short.q \
                -l h_rt=48:00:00  \
                -l mfree=10G \
                -l gpfsstate=0 \
                -pe serial {threads} \
                -V -cwd \
                -o $LOGDIR \
                -e $LOGDIR \
                -S /bin/bash" \
        --jobs $JOBNUM \
        --latency-wait $WAITTIME \
        --restart-times $RETRY \
        -p \
        -k $1 $2
