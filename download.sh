#!/bin/bash

echo "first $1";
echo "last $2";
server="uqyche38@goliath.labs.eait.uq.edu.au"
server_prefix="~/macondo/0.5r_10.0Ga_0.3gap/"
local_prefix="/home/pzhang/chen/macondo/0.5r_10.0Ga_0.3gap/"
file_prefix="test_mvbed_c_"
echo $server_line
for i in $(seq $1 $2);
do
	scp -r $server:$server_prefix$file_prefix$(printf %04d $i).h5 $local_prefix
	# scp -r uqyche38@goliath.labs.eait.uq.edu.au:~/macondo/0.5r_10.0Ga_0.3gap/test_mvbed_c_$(printf %04d $i).xmf ;
	scp -r $server:$server_prefix$file_prefix$(printf %04d $i).xmf $local_prefix
	python process.py $local_prefix$file_prefix $i $i
done;
