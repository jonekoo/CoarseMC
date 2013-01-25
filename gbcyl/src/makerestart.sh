#!/bin/sh

begin=$1
end=$2
i=$begin

while [ $i -lt $end ]; 
do
  test -f inputconfiguration.$i && \
  mv inputconfiguration.$i inputconfiguration.$i~
  cp restartconfiguration.$i inputconfiguration.$i

  test -f inputparameters.$i && \
  mv inputparameters.$i inputparameters.$i~
  cp restartparameters.$i inputparameters.$i

  test -f mtstate.$i && \
  mv mtstate.$i mtstate.$i~
  cp restartmtstate.$i mtstate.$i

  i=`expr $i + 1`
done