#!/usr/bin/bash
#
COORD="${1:-sph}"
PROB="${2:-blast}"
echo "Building with coord: "
echo "$COORD"
echo " and "
echo "prob: "
echo $PROB
echo "\n"
cd $AD
if [[ "$COORD" == "cart" ]]; then
  python3 configure.py --prob=$PROB --coord=cartesian
fi
if [[ "$COORD" == "sph" ]]; then
    python3 configure.py --prob=$PROB --coord=spherical_polar
fi
make clean
make
