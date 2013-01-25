#!/bin/sh
#
# Runs the script given as argument to option -s over arguments of this script.
# The arguments should be file names. 
#
# Author: 
# Jouni Karjalainen 
# Department of Physical Sciences
# University of Oulu 
#
# Much of this script owes to the scripts presented by Ken Steube 
# (http://www.uq.edu.au/~uqksteub/) in his UNIX Bourne Shell Scripting 
# tutorial which used to be available in the WWW. 
#
#
# Print usage: loop-script-over-args.sh -h
#

print_usage()
{
  echo "Usage: loop-script-over-args.sh [-s script || -h] [arg1 arg2 ...]"
}

# Initialize our variables so we don't inherit values
# from the environment

opt_s=''
opt_h=''

# Parse the command-line options
while getopts 's:h' option
do
case "$option" in
    "s")opt_s="$OPTARG"
;;
    "h")opt_h="1"
;;
    ?)print_usage
exit 1
;;
esac
done

shift `expr $OPTIND - 1`

if [ "$opt_s" != "" ]
then
  if [ "$*" != "" ]
  then
  for arg in "$@"
  do
      sh "$opt_s" "$arg"
  done
  fi
else
  if [ "$opt_h" != "" ]
  then
    print_usage
  fi
fi
