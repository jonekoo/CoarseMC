#!/usr/bin/awk -f
#
# For each column in a file, concatenates it to a string and feeds it to the 
# given command. The concatenation turns the column into a row. So for example
# the if you want to turn your columns into rows (matrix transponation)  just 
# use -v com=cat. 
# 
# Author:
# Jouni Karjalainen 
# Department of Physical Sciences
# University of Oulu
#
# Usage: awk -v com=command -f pipe-columns.awk filename
#
{
  for(i=1; i<=NF; i=i+1) 
  {
    result[i] = result[i] " " $i
  }
}
END {
  for(i=1; i<=NF; i=i+1) 
  {
    system("echo " result[i] "|" com)
  }
}

