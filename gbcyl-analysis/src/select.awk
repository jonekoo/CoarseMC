#!/usr/bin/awk -f
#
# Prints the configurations i to which start <= i < stop
#
# Usage example: 
# cat simdata.out| select.awk -v start=10 -v stop=11
# Prints the 10th configuration. 
#
BEGIN {i=0}
{
  if ($0 ~ /R/) i=i+1;
  if (i>=start && i<stop) print $0
}