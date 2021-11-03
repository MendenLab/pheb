#!/bin/bash

cancertype=$1

elements=$(ls metadata/results_dmp_new/$cancertype/$cancertype.meth.*.csv)
elements=($elements)

if [ "$2" = "test" ]
then
	elementss=${elements[@]:0:10}
fi 

if [ "$2" != "test" ]
then
        elementss=${elements[@]}
fi

elementss=($elementss)
echo Performing for following files...
echo ${elementss[@]}


for element in ${elementss[@]}
do
	echo $element
	combined-pvalues-master/cpv/comb-p pipeline -c 5 --seed 0.001 --dist 200 -p $element.combp --region-filter-p 0.05 --anno mm9 $element
done

