#!/bin/bash

echo "first $1";
echo "last $2";

for i in $(seq $1 $2);
do
	scp -r uqyche38@goliath.labs.eait.uq.edu.au:~/macondo/test_mvbed_c_$(printf %04d $i).h5 ~/chen/macondo;
	scp -r uqyche38@goliath.labs.eait.uq.edu.au:~/macondo/test_mvbed_c_$(printf %04d $i).xmf ~/chen/macondo;
done;
