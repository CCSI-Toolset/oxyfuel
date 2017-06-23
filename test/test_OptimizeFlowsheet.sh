#!/bin/sh

THIS_SCRIPT=`basename $0`
SCRIPT_DIR=`pwd`
LOG=${SCRIPT_DIR}/${THIS_SCRIPT}.log

WORK_DIR=`cd ${SCRIPT_DIR}/../worker0; pwd`
OUTPUT=${WORK_DIR}/OptimizeFlowsheet.lst

EXPECTED_VAL=0.2853

die () {
    echo "${THIS_SCRIPT}: $1";
    exit "$2";
}

cd $WORK_DIR || die "no worker0 dir" 1

gams ../OptimizeFlowsheet.gms LogOption=3 > $LOG
# Will return non-zero (7) due to lack of license even if run
# succeeds, so no bother checking status

# Check for output file
[ -e ${OUTPUT} ] || die "output file (${OUTPUT}) not found" 1

# Check results
FOUND_VAL=`grep 'VAR Z9' $OUTPUT | tail -1 | awk '{print $5}'`
if [ "${FOUND_VAL}" = "${EXPECTED_VAL}" ]; then
    exit 0;
else
    msg="Found ${FOUND_VAL}, expected ${EXPECTED_VAL}\n
For details see one or more of the following:\n
  $OUTPUT\n
  $LOG";
    die "${msg}" 1
fi
