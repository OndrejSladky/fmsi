#!/usr/bin/env bash
BIN=bin
PROG=../fmsi
TESTS=testfiles

mkdir -p $BIN

# The commands of FMSI 0.4, on its index of several files.
LEGACY="$PROG legacy"

$LEGACY index $TESTS/integration_a.fa
$LEGACY index $TESTS/integration_b.fa

$LEGACY merge -p $TESTS/integration_a.fa -p $TESTS/integration_b.fa -r $BIN/merged.fa

$LEGACY query -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a.txt 2> /dev/null
$LEGACY lookup -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a_hash.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_a.fa > $BIN/a_xor.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt $TESTS/integration_b.fa > $BIN/b.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_b.fa > $BIN/b_xor.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt $BIN/merged.fa > $BIN/merged.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $BIN/merged.fa > $BIN/merged_xor.txt 2> /dev/null

#$LEGACY normalize -k 3 -p $BIN/merged.fa -s -l > $BIN/merged_normalized.fa 2> /dev/null
$LEGACY normalize -k 3 -s $BIN/merged.fa > $BIN/merged_normalized2.fa 2> /dev/null



diff $TESTS/result_a_complements.txt $BIN/a.txt || exit 1
echo "a.txt OK"
diff $TESTS/result_a_complements_hash.txt $BIN/a_hash.txt || exit 1
echo "a_hash.txt OK"
diff $TESTS/result_b_complements.txt $BIN/b.txt || exit 1
echo "b.txt OK"
diff $TESTS/result_b_complements_xor.txt $BIN/b_xor.txt || exit 1
echo "b_xor.txt OK"
diff $TESTS/result_merged_complements.txt $BIN/merged.txt || exit 1
echo "merged.txt OK"

#diff $TESTS/result_normalized.txt $BIN/merged_normalized.fa
diff $TESTS/result_normalized2.txt $BIN/merged_normalized2.fa || exit 1
echo "merged_normalized2.fa OK"

# FMSI 2: one index file, with the same answers as FMSI 0.4, in both presets and every layout.
for preset in extra-small small fast; do
  $PROG index --$preset -l -k 3 -o $BIN/a_$preset $TESTS/integration_a.fa 2> /dev/null || exit 1
  [ -f $BIN/a_$preset.fmsi ] && [ -f $BIN/a_$preset.fmsi.klcp ] || exit 1
  $PROG info --check $BIN/a_$preset.fmsi > $BIN/a_$preset.info || exit 1
  grep -q "^preset	$preset$" $BIN/a_$preset.info || exit 1
  grep -q "^k	3$" $BIN/a_$preset.info || exit 1
  for layout in "" "-c 64 -m plain" "-c 128 -m packed" "-c 256 -m rrr"; do
    for streaming in "" "-S" "-S -L"; do
      $PROG query $layout $streaming -q $TESTS/queries.txt $BIN/a_$preset.fmsi > $BIN/a2.txt 2> /dev/null || exit 1
      diff $TESTS/result_a_complements.txt $BIN/a2.txt || exit 1
    done
  done
  # The index is also found by its name without the extension.
  $PROG query -q $TESTS/queries.txt $BIN/a_$preset > $BIN/a2.txt 2> /dev/null || exit 1
  diff $TESTS/result_a_complements.txt $BIN/a2.txt || exit 1
done
echo "FMSI 2 membership OK"

# Dictionary queries: the hashes of the --small index are those of FMSI 0.4.
$PROG lookup --extra-small -q $TESTS/queries.txt $BIN/a_small.fmsi > $BIN/a2_hash.txt 2> /dev/null || exit 1
diff $TESTS/result_a_complements_hash.txt $BIN/a2_hash.txt || exit 1
# On a min-one mask, lookup is a minimal perfect hash and access inverts it.
MINONE=data/GCF_test2_ones.fna.ms.k31
$PROG index -o $BIN/minone $MINONE 2> /dev/null || exit 1
$PROG info $BIN/minone.fmsi > $BIN/minone.info || exit 1
grep -q "^mask_type	min-one" $BIN/minone.info || exit 1
N=$(grep "^mask_ones" $BIN/minone.info | cut -f2)
seq 0 $((N - 1)) > $BIN/hashes.txt
$PROG access -q $BIN/hashes.txt $BIN/minone.fmsi > $BIN/kmers.txt 2> /dev/null || exit 1
awk '{ print ">" NR; print $0 }' $BIN/kmers.txt > $BIN/kmers.fa
$PROG lookup -q $BIN/kmers.fa $BIN/minone.fmsi 2> /dev/null | cut -f2 > $BIN/hashes_back.txt || exit 1
diff $BIN/hashes.txt $BIN/hashes_back.txt || exit 1
# The type of a mask that sets a k-mer in both orientations is found out on both strands only.
OTHER=data/GCF_009858895.2_ASM985889v3_genomic.fna.ms.k13
$PROG index --fast -o $BIN/other $OTHER 2> /dev/null || exit 1
$PROG info $BIN/other.fmsi | grep -q "^mask_type	other$" || exit 1
echo "FMSI 2 dictionary OK"

# The mask stored along the superstring (--compact-file): the same answers; export gives the superstring back.
for preset in small fast; do
  $PROG index --$preset -z -o $BIN/compact_$preset $MINONE 2> /dev/null || exit 1
  $PROG info --check $BIN/compact_$preset.fmsi | grep -q "^mask_order	text$" || exit 1
  $PROG index --$preset -o $BIN/plain_$preset $MINONE 2> /dev/null || exit 1
  for layout in "" "--extra-small"; do
    $PROG lookup $layout -S -q $BIN/kmers.fa $BIN/compact_$preset.fmsi > $BIN/compact.txt 2> /dev/null || exit 1
    $PROG lookup $layout -S -q $BIN/kmers.fa $BIN/plain_$preset.fmsi > $BIN/plain.txt 2> /dev/null || exit 1
    diff $BIN/plain.txt $BIN/compact.txt || exit 1
  done
  $PROG export $BIN/compact_$preset.fmsi > $BIN/exported.fa 2> /dev/null || exit 1
  diff <(sed -n 2p $MINONE) <(sed -n 2p $BIN/exported.fa) || exit 1
done
[ "$(wc -c < $BIN/compact_small.fmsi)" -lt "$(wc -c < $BIN/plain_small.fmsi)" ] || exit 1
echo "FMSI 2 compact file and export OK"

echo "All tests passed"

