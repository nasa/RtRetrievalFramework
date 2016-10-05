#!/bin/sh
#sounding_id            lat         lon   surf_type        gain        Site   xco2_true
#YYYYMMDDhhmmss         deg         deg                    none        none         ppm
output_file=quick_look_comparison.txt
tmp_file=`mktemp`
tail -n 10000 valset_model_10k_Enhanced_Regression_20130327.txt > $tmp_file
tail -n 2600 15000-TCCON-nearest-match-to-ACOS-sounding_FULL_Glint.csv | sed 's|Park Falls|Park_Falls|' | sed 's|,|\t|g' | sed 's|"||g' | awk '{print $1"\t"$2"\t"$3"\tocean\tH\t"$4"\t"$5}' >> $tmp_file
tail -n 10400 15000-TCCON-nearest-match-to-ACOS-sounding_FULL_Hgain.csv | sed 's|Park Falls|Park_Falls|' | sed 's|,|\t|g' | sed 's|"||g' | awk '{print $1"\t"$2"\t"$3"\tland\tH\t"$4"\t"$5}' >> $tmp_file
tail -n 2000 15000-TCCON-nearest-match-to-ACOS-sounding_FULL_Mgain.csv  | sed 's|Park Falls|Park_Falls|' | sed 's|,|\t|g' | sed 's|"||g' | awk '{print $1"\t"$2"\t"$3"\tland\tM\t"$4"\t"$5}' >> $tmp_file

head -n 2 valset_model_10k_Enhanced_Regression_20130327.txt > $output_file
sort $tmp_file >> $output_file

rm $tmp_file
