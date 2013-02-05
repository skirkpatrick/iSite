#!/bin/bash

mkdir results

for i in $(ls | grep results.)
do
    cp average.rb $i
    cp average_dist.rb $i
    cd $i
    ./average.rb result.* $i
    ./average_dist.rb dist* dist_$i
    mv $i ../results
    mv dist_$i ../results
    rm average.rb
    rm average_dist.rb
    cd ..
done

cp collect.rb results

cd results

./collect.rb results.* result
rm collect.rb
