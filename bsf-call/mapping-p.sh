#!/bin/sh

date_format="%Y-%m-%d %H:%M:%S"

echo `date +"$date_format"` $0 start.

read_dir=$1
result_dir=$2
refgenome=$3
map_probs_opt=$4
read_name1=$5
read_format1=$6
lastal_opt1=$7
read_name2=$8
read_format2=$9
lastal_opt2=${10}

echo Host: `hostname`
echo Arguments:
echo "  Read directory:              $read_dir"
echo "  Result directory:            $result_dir"
echo "  Reference genome:            $refgenome"
echo "  Options for last-pair-probs: $map_probs_opt"
echo "  Read name          (read1):  $read_name1"
echo "  Read format        (read1):  $read_format1"
echo "  Options for lastal (read1):  $lastal_opt1"
echo "  Read name          (read2):  $read_name2"
echo "  Read format        (read2):  $read_format2"
echo "  Options for lastal (read2):  $lastal_opt2"
echo

echo `date +"$date_format"` gzip read files start. 1>&2

pushd $read_dir
echo `date +"$date_format"` gzip "$read_name1"_1."$read_format1"
gzip "$read_name1"_1."$read_format1" &
pid1=$!

echo `date +"$date_format"` gzip "$read_name2"_2."$read_format2"
gzip "$read_name2"_2."$read_format2" &
pid2=$!
popd

wait $pid1
wait $pid2

echo `date +"$date_format"` gzip read files done. 1>&2
echo `date +"$date_format"` gzip read files done.

echo `date +"$date_format"` 'lastal (1st read, forward) start.' 1>&2
echo `date +"$date_format"` zcat $read_dir/"$read_name1"_1."$read_format1".gz '|' lastal -p bisulfite_f.mat -s1 $lastal_opt1 $refgenome.f '- >' $result_dir/"$read_name1"_1_f.maf
zcat $read_dir/"$read_name1"_1."$read_format1".gz | lastal -p bisulfite_f.mat -s1 $lastal_opt1 $refgenome.f - > $result_dir/"$read_name1"_1_f.maf
echo `date +"$date_format"` 'lastal (1st read, forward) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'lastal (1st read, reverse) start.' 1>&2
echo `date +"$date_format"` zcat $read_dir/"$read_name1"_1."$read_format1".gz '|' lastal -p bisulfite_r.mat -s0 $lastal_opt1 $refgenome.r '- >' $result_dir/"$read_name1"_1_r.maf
zcat $read_dir/"$read_name1"_1."$read_format1".gz | lastal -p bisulfite_r.mat -s0 $lastal_opt1 $refgenome.r - > $result_dir/"$read_name1"_1_r.maf
echo `date +"$date_format"` 'lastal (1st read, reverse) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'last-merge-batches (1st read) start.' 1>&2
echo `date +"$date_format"` last-merge-batches $result_dir/"$read_name1"_1_f.maf $result_dir/"$read_name1"_1_r.maf '>' $result_dir/"$read_name1"_1.maf
last-merge-batches $result_dir/"$read_name1"_1_f.maf $result_dir/"$read_name1"_1_r.maf > $result_dir/"$read_name1"_1.maf
echo `date +"$date_format"` 'last-merge-batches (1st read) done.' exit-status: $? 1>&2

echo `date +"$date_format"` delete temporary files start. 1>&2
rm $result_dir/"$read_name1"_1_f.maf $result_dir/"$read_name1"_1_r.maf
echo `date +"$date_format"` delete temporary files done. exit-status: $? 1>&2

echo `date +"$date_format"` 'lastal (2nd read, forward) start.' 1>&2
echo `date +"$date_format"` zcat $read_dir/"$read_name2"_2."$read_format2".gz '|' lastal -p bisulfite_f.mat -s1 $lastal_opt2 $refgenome.f '- >' $result_dir/"$read_name2"_2_f.maf
zcat $read_dir/"$read_name2"_2."$read_format2".gz | lastal -p bisulfite_f.mat -s1 $lastal_opt2 $refgenome.f - > $result_dir/"$read_name2"_2_f.maf
echo `date +"$date_format"` 'lastal (2nd read, forward) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'lastal (2nd read, reverse) start.' 1>&2
echo `date +"$date_format"` zcat $read_dir/"$read_name2"_2."$read_format2".gz '|' lastal -p bisulfite_r.mat -s0 $lastal_opt2 $refgenome.r '- >' $result_dir/"$read_name2"_2_r.maf
zcat $read_dir/"$read_name2"_2."$read_format2".gz | lastal -p bisulfite_r.mat -s0 $lastal_opt2 $refgenome.r - > $result_dir/"$read_name2"_2_r.maf
echo `date +"$date_format"` 'lastal (2nd read, reverse) done.' exit-status: $? 1>&2

echo `date +"$date_format"` 'last-merge-batches (2nd read) start.' 1>&2
echo `date +"$date_format"` last-merge-batches $result_dir/"$read_name2"_2_f.maf $result_dir/"$read_name2"_2_r.maf '>' $result_dir/"$read_name2"_2.maf
last-merge-batches $result_dir/"$read_name2"_2_f.maf $result_dir/"$read_name2"_2_r.maf > $result_dir/"$read_name2"_2.maf
echo `date +"$date_format"` 'last-merge-batches (2nd read) done.' exit-status: $? 1>&2

echo `date +"$date_format"` delete temporary files start. 1>&2
rm $result_dir/"$read_name2"_2_f.maf $result_dir/"$read_name2"_2_r.maf
echo `date +"$date_format"` delete temporary files done. exit-status: $? 1>&2

echo `date +"$date_format"` last-pair-probs start. 1>&2
echo `date +"$date_format"` last-pair-probs $map_probs_opt $result_dir/"$read_name1"_1.maf $result_dir/"$read_name2"_2.maf '>' $result_dir/"$read_name1".maf
last-pair-probs $map_probs_opt $result_dir/"$read_name1"_1.maf $result_dir/"$read_name2"_2.maf > $result_dir/"$read_name1".maf
echo `date +"$date_format"` last-pair-probs done. exit-status: $? 1>&2

echo `date +"$date_format"` delete temporary files start. 1>&2
rm $result_dir/"$read_name1"_1.maf $result_dir/"$read_name2"_2.maf
echo `date +"$date_format"` delete temporary files done. exit-status: $? 1>&2

echo `date +"$date_format"` 'gzip result file ('$result_dir/"$read_name1".maf') start.' 1>&2
gzip $result_dir/"$read_name1".maf
echo `date +"$date_format"` gzip result file done. exit-status: $? 1>&2

echo `date +"$date_format"` $0 done.

exit 0
