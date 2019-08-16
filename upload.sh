#!/bin/bash

echo "first $1";
echo "last $2";
server="uqyche38@goliath.labs.eait.uq.edu.au"
#server="uqyche38@awoonga.qriscloud.org.au"
server_prefix="~/macondo/5e3Re_27.0Ga_0.3gap_a/"
#server_prefix="~/macondo/small/1/"
#server_prefix="/30days/uqyche38/0.05r_10.0Ga_0.3gap/"
IFS=$'\n'
#local_prefix=/media/pzhang/Elements/move-bed-tmp/macondo/0.05r_10.0Ga_0.3gap/
#local_prefix=/media/pzhang/My\ Book/dune_shape/small/1/
local_prefix=/PZ_Q/move-bed-tmp/5e3Re_27.0Ga_0.3gap_a/
file_prefix="test_mvbed_c_"
echo $server_line
for i in $(seq $1 $2);
do
	scp -r $local_prefix$file_prefix$(printf %04d $i).h5 $server:$server_prefix
	scp -r $local_prefix$file_prefix$(printf %04d $i).xmf $server:$server_prefix
	
done;
