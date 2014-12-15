#!/bin/sh

date_format="%Y-%m-%d %H:%M:%S"

echo `date +"$date_format"` $0 start.

read_dir=$1
result_dir=$2
refgenome=$3
map_probs_opt=$4
read_name=$5
read_format=$6
lastal_opt=$7

echo Host: `hostname`
echo Arguments:
echo "  Read directory:             $read_dir"
echo "  Result directory:           $result_dir"
echo "  Reference genome:           $refgenome"
echo "  Options for last-map-probs: $map_probs_opt"
echo "  Read name:                  $read_name"
echo "  Read format:                $read_format"
echo "  Options for lastal:         $lastal_opt"
echo

echo `date +"$date_format"` gzip read file start. 1>&2
pushd $read_dir
echo `date +"$date_format"` gzip "$read_name"_1."$read_format"; gzip "$read_name"_1."$read_format"
popd
echo `date +"$date_format"` gzip read file done. 1>&2


echo `date +"$date_format"` 'lastal (forward) start.' 1>&2
echo `date +"$date_format"` zcat "$read_dir/$read_name"_1."$read_format".gz '|' lastal -p bisulfite_f.mat -s1 $lastal_opt "$refgenome".f '- >' "$result_dir/$read_name"_1_f.maf
zcat "$read_dir/$read_name"_1."$read_format".gz | lastal -p bisulfite_f.mat -s1 $lastal_opt "$refgenome".f - > "$result_dir/$read_name"_1_f.maf
echo `date +"$date_format"` 'lastal (forward) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'lastal (reverse) start.' 1>&2
echo `date +"$date_format"` zcat "$read_dir/$read_name"_1."$read_format".gz '|' lastal -p bisulfite_r.mat -s0 $lastal_opt "$refgenome".r '- >' "$result_dir/$read_name"_1_r.maf
zcat "$read_dir/$read_name"_1."$read_format".gz | lastal -p bisulfite_r.mat -s0 $lastal_opt "$refgenome".r - > "$result_dir/$read_name"_1_r.maf
echo `date +"$date_format"` 'lastal (reverse) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'last-merge-batches start.' 1>&2
echo `date +"$date_format"` last-merge-batches "$result_dir/$read_name"_1_f.maf "$result_dir/$read_name"_1_r.maf '>' "$result_dir/$read_name"_1.maf
last-merge-batches "$result_dir/$read_name"_1_f.maf "$result_dir/$read_name"_1_r.maf > "$result_dir/$read_name"_1.maf
echo `date +"$date_format"` 'last-merge-batches done.' exit-status: $? 1>&2

echo `date +"$date_format"` delete temporary files start. 1>&2
rm "$result_dir/$read_name"_1_f.maf "$result_dir/$read_name"_1_r.maf
echo `date +"$date_format"` delete temporary files done. exit-status: $? 1>&2

echo `date +"$date_format"` 'last-map-probs start.' 1>&2
echo `date +"$date_format"` last-map-probs $map_probs_opt "$result_dir/$read_name"_1.maf '>' "$result_dir/$read_name".maf
last-map-probs $map_probs_opt "$result_dir/$read_name"_1.maf > "$result_dir/$read_name".maf
echo `date +"$date_format"` 'last-map-probs done.' exit-status: $? 1>&2

echo `date +"$date_format"` delete temporary files start. 1>&2
rm "$result_dir/$read_name"_1.maf
echo `date +"$date_format"` delete temporary files done. exit-status: $? 1>&2

echo `date +"$date_format"` 'gzip result file ('$result_dir/"$read_name".maf') start.' 1>&2
gzip "$result_dir/$read_name".maf
echo `date +"$date_format"` gzip result file done. exit-status: $? 1>&2

echo `date +"$date_format"` $0 done.
