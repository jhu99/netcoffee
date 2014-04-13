#!/bin/bash
for ((i=0;i<4;i++))
do
	mv net$i.tab ppi$i.tab
	for ((j=i;j<4;j++))
	do
		mv net$i-net$j.cf ppi$i-ppi$j.cf
	done
done
