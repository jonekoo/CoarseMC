#!/usr/bin/awk -f
#
# Usage example: 
# cat simdata.out| select.awk -v pattern=R -v start=10 -v stop=11
# Prints the 10th configuration. 
#
BEGIN {i=0}
{
  if ($0 ~ /R/) i=i+1;
  if (i>=start && i<stop) print $0
}