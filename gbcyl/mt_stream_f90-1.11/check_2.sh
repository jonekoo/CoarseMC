#!/bin/sh

echo "check_stream_2 : consistency check #2."
echo "check_stream_2 : checking multiple stream consistency with huge jump."
 
#tempfile=check_stream_2_`date "+%Y-%m-%dT%H:%M:%S"`.out
#./check_stream_2 > $tempfile


diff test_2.out sample_2.out
if [ $? -ne 0 ]; then
  echo "check_stream_2 : Check Error."
  exit -1
else
  echo "check_stream_2 : Check OK."
  exit 0
fi
#rm tempfile
