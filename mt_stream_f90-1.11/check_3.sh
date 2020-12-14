#!/bin/sh

echo "check_stream_3 : consistency check #3."
echo "check_stream_3 : checking state save/read."
./check_stream_3 > test_3.out
grep '^S' test_3.out | sed 's/S//g' > save.dat
grep '^L' test_3.out | sed 's/L//g' > read.dat
diff save.dat read.dat

if [ $? -ne 0 ]; then
  echo "check_stream_3 : Check Error."
  exit -1
else
  echo "check_stream_3 : Check OK."
  exit 0
fi

