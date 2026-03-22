#!/bin/bash
# Optimized 90d factoring script
# Uses YAFU GNFS with GGNFS sievers
# Attempts to factor a given 90-digit semiprime in under 300s single-core
#
# Usage: timeout 295 bash library/factor_90d.sh <N>

N=$1
if [ -z "$N" ]; then
    echo "Usage: $0 <number>"
    exit 1
fi

YAFU=/tmp/agent-factoring-4/yafu/yafu
GGNFS_DIR=/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/

WORKDIR=$(mktemp -d /tmp/factor90_XXXXXX)
cd $WORKDIR

cat > yafu.ini << EOF
ggnfs_dir=$GGNFS_DIR
EOF

echo "nfs($N)" | $YAFU -threads 1 -seed 42 -xover 85 2>&1

EXIT=$?
rm -rf $WORKDIR
exit $EXIT
