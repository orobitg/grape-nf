#!/bin/bash

set -e
set -u

IN=${1}
LEN=${2}
OUTFILE=${3}

BAMNAME=$(basename ${OUTFILE})
OUTDIR=$(dirname ${OUTFILE})


if [ ${LEN} -ge 2 ]; then

	$HOME/CRG/Projects/Selenoprofiles/grape-nf/bin/MergeSamFiles.jar ${IN} O=${BAMNAME}.bam USE_THREADING=true
else 
	mv ${IN} ${BAMNAME}.bam
fi

#	BAM=${BAMNAME}.bam
#	rm -f ${OUTDIR}/${BAM}
#	mv ${BAM} ${OUTDIR}/${BAM}
#	ln -s ${OUTDIR}/${BAM}

#elif [ $l -eq 1 ]; then
#	BAM=`ls ${OUTDIR} | grep '.bam' | head -$l | tail -1`
#	cp ${OUTDIR}/${BAM} ${BAM}
#fi
