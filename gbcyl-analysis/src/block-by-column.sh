#!/bin/sh
#
# Author: 
# Jouni Karjalainen 
# Department of Physical Sciences
# University of Oulu 
#
# The command line option parsing of this script owes to the scripts presented
# by Ken Steube (http://www.uq.edu.au/~uqksteub/) in his UNIX Bourne Shell 
# Scripting tutorial which used to be available in the WWW. 
#
#
# Print usage: block-by-column.sh -h
#

print_usage()
{
  echo "Usage: block-by-column.sh [-f file || -h] [arg1 arg2 ...]"
}

# Initialize our variables so we don't inherit values
# from the environment

opt_f=''
opt_h=''

# Parse the command-line options
while getopts 'f:h' option
do
case "$option" in
    "f")opt_f="$OPTARG"
;;
    "h")opt_h="1"
;;
    ?)print_usage
exit 1
;;
esac
done

shift `expr $OPTIND - 1`

if [ "$opt_f" != "" ]
then
  i="1"
  ncolumns=`awk '{if(1==NR) {print NF}}' $opt_f`
  while [ $i -le $ncolumns ]
  do 
    cat $opt_f | awk -v column=$i '{print $column}' | blocked_error
    i=`expr $i + 1`
  done 
else
  if [ "$opt_h" != "" ]
  then
    print_usage
  fi
fi

