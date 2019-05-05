#!/bin/bash

server_prefix="/30days/uqyche38/0.5r_10.0Ga_0.3gap/"
file_prefix="test_mvbed_c_cti1_"
echo "First $1";
echo "Second $2";
for i in $(seq $1 $2);
do
	rm $server_prefix$file_prefix$(printf %04d $i).h5
	rm $server_prefix$file_prefix$(printf %04d $i).xmf
done
