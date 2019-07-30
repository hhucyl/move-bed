#!/bin/bash

echo "first $1";
echo "last $2";
#server="uqyche38@goliath.labs.eait.uq.edu.au"
server="uqyche38@awoonga.qriscloud.org.au"
#server_prefix="~/macondo/0.2r_10.0Ga_0.3gap/"
#server_prefix="~/macondo/small/1/"
server_prefix="/30days/uqyche38/5e3Re_26.2Ga_0.3gap_a/"
IFS=$'\n'
#local_prefix=/media/pzhang/Elements/move-bed-tmp/macondo/0.002r_20.0Ga_0.3gap/
local_prefix=/media/pzhang/Elements/move-bed-tmp/macondo/5e3Re_26.2Ga_0.3gap_a/
#local_prefix=/media/pzhang/My\ Book/dune_shape/small/1/
file_prefix="test_mvbed_c_"
echo $server_line
for i in $(seq $1 $2);
do
	scp -r $server:$server_prefix$file_prefix$(printf %04d $i).h5 $local_prefix
	# scp -r uqyche38@goliath.labs.eait.uq.edu.au:~/macondo/0.5r_10.0Ga_0.3gap/test_mvbed_c_$(printf %04d $i).xmf ;
	scp -r $server:$server_prefix$file_prefix$(printf %04d $i).xmf $local_prefix
	python process.py $local_prefix$file_prefix $i $i
done;
