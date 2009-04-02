#!/usr/bin/awk -f
#
# Usage example: 
# cat simdata.out| periodic_select.awk -v period=10
# Prints every 10th configuration. 
#
BEGIN {i=0}
{
  if ($0 ~ /R/) i=i+1;
  if (i % period == 0) print $0
}