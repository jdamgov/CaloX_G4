#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters.  Run: runLocalSim.sh pi+ 30"
    exit 1
fi

cd /home/jdamgov/CaloX/sim

PE=`echo  "$1$2"  | sed  "s/\-$2/M$2/g" |sed  "s/\+$2/P$2/g"`

#rm -rf $PE
mkdir -p $PE 
cp paramBatch_tmp.mac $PE/paramBatch.mac
sed -i "s/XXXXX/$1/g" $PE/paramBatch.mac
sed -i "s/YYYYY/$2/g" $PE/paramBatch.mac

cp runBatch_tmp.sh $PE/runBatch_$3.sh
sed -i "s/XXXXX/$1/g" $PE/runBatch_$3.sh
sed -i "s/YYYYY/$2/g" $PE/runBatch_$3.sh
sed -i "s/ZZZZZ/$3/g" $PE/runBatch_$3.sh

sed -i "s/\-$2/M$2/g" $PE/runBatch_$3.sh
sed -i "s/\+$2/P$2/g" $PE/runBatch_$3.sh

sleep 10
$PE/runBatch_$3.sh > $PE/runBatch_$3.out &
