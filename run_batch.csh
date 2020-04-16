#!/bin/bash
TEMPFILE=./cnt.tmp
#echo 0 > $TEMPFILE


while [ Your != "done" ]
do
  export nums=`ps aux |grep exampleB4b |wc -l`
  echo "nums" $nums
  if [ $nums != "8" ] && [ $nums != "9" ] && [ $nums != "7" ]; then
      COUNTER=$[$(cat $TEMPFILE) + 1]
      #./runLocalSimL.sh e+ 150 $COUNTER
      ./runLocalSimL.sh pi+ 30 $COUNTER
      echo $COUNTER > $TEMPFILE
      echo $COUNTER
  fi
  sleep 10
done
 
