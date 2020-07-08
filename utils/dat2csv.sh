#!/bin/bash

# This script with prepare *.csv file with columns as follows
# for a range of bands [nbmin..nbmax]
#
# output file *.csv
# k [rad/bohr] along the path, E(nbmin), E(nbmin+1), ..., E(nbmax), minv(nbmin), minv(nbmin+1), ..., minv(nbmax)
# ...
# here E are energy eigenvalues (eV) relative to the Fermi energy set diring band structure calculation
# minv is the inverse effective mass component selected below
#
# (c) Oleg Rubel (May 2020)

nbmin=21
nbmax=22
case=${PWD##*/}
#fname='minv_c-up'; col=2; suffix='' # conductivity eff. mass; col = column in the *dat file
#fname='minv_d'; col=2; suffix='' # DOS eff. mass; col = column in the *dat file
fname='minv_ij-up'; col=2; suffix='_xx' # eff. mass xx
#fname='minv_ij-up'; col=3; suffix='_yy' # eff. mass yy
#fname='minv_ij'; col=4; suffix='_zz' # eff. mass zz



# get column of dk on the path

rm tmp1
sed -n -E "/^  bandindex:[[:space:]]+1$/,/^  bandindex:[[:space:]]+2$/p" ${case}.spaghetti_ene | sed '/^  bandindex/d' | awk {'print $4'} > tmp1 # [[:space:]]+ will match arbitrary number of spaces


# get band energies [eV]

for i in $(seq $nbmin $nbmax) # loop over bands in the range
do
  echo "band $i"
  ip1=$((i + 1))
  echo "ip1 = $ip1"
  sed -n -E "/^  bandindex:[[:space:]]+${i}$/,/^  bandindex:[[:space:]]+${ip1}$/p" ${case}.spaghetti_ene | sed '/^  bandindex/d' | awk {'print $5'} > tmp2
  paste -d, tmp1 tmp2 > tmp3 # past in CSV format
  mv tmp3 tmp1
done

# get columns of effective masses

for i in $(seq $nbmin $nbmax) # loop over bands in the range
do
  echo "band $i"
  grep -E "^ +${i} " ${fname}.dat | awk {'print $'$col} > tmp2 # here ' +' means arbitrary number of spaces
  paste -d, tmp1 tmp2 > tmp3 # past in CSV format
  mv tmp3 tmp1
done
rm tmp2 # leave tmp1

# prepare CSV file

echo "# k [rad/bohr],energy [eV] for bands ${nbmin}-${nbmax}, (m0/m*) for bands ${nbmin}-${nbmax}" > ${fname}${suffix}-b${nbmin}-${nbmax}.csv
cat tmp1 >> ${fname}${suffix}-b${nbmin}-${nbmax}.csv
rm tmp1
echo "Results are stored in file ${fname}${suffix}-b${nbmin}-${nbmax}.csv"

# Display tick label that are usefull for gnuplot

echo "Here are ticklabel that should be used in gnuplot under 'set xtics ...'"
grep -E "ticklabel +[0-9]+ " -B 1 *.bands.agr
