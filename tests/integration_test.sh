#!/usr/bin/env bash
BIN=bin
PROG=../fmsi
TESTS=testfiles

mkdir -p $BIN

$PROG index $TESTS/integration_a.fa
$PROG index $TESTS/integration_b.fa

$PROG merge -p $TESTS/integration_a.fa -p $TESTS/integration_b.fa -r $BIN/merged.fa

$PROG query -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a.txt 2> /dev/null
$PROG lookup -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a_hash.txt 2> /dev/null
$PROG query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_a.fa > $BIN/a_xor.txt 2> /dev/null
$PROG query -k 3 -q $TESTS/queries.txt $TESTS/integration_b.fa > $BIN/b.txt 2> /dev/null
$PROG query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_b.fa > $BIN/b_xor.txt 2> /dev/null
$PROG query -k 3 -q $TESTS/queries.txt $BIN/merged.fa > $BIN/merged.txt 2> /dev/null
$PROG query -k 3 -q $TESTS/queries.txt -f xor $BIN/merged.fa > $BIN/merged_xor.txt 2> /dev/null

#$PROG normalize -k 3 -p $BIN/merged.fa -s -l > $BIN/merged_normalized.fa 2> /dev/null
$PROG normalize -k 3 -s $BIN/merged.fa > $BIN/merged_normalized2.fa 2> /dev/null



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

echo "All tests passed"

