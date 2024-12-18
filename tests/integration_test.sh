#!/usr/bin/env bash
BIN=bin
PROG=../fmsi
TESTS=testfiles

mkdir -p $BIN

$PROG index -p $TESTS/integration_a.fa
$PROG index -p $TESTS/integration_b.fa

$PROG merge -p $TESTS/integration_a.fa -p $TESTS/integration_b.fa -r $BIN/merged.fa

$PROG query -F -k 3 -p $TESTS/integration_a.fa -q $TESTS/queries.txt > $BIN/a.txt 2> /dev/null
$PROG query -F -k 3 -p $TESTS/integration_a.fa -q $TESTS/queries.txt -f xor > $BIN/a_xor.txt 2> /dev/null
$PROG query -F -k 3 -p $TESTS/integration_b.fa -q $TESTS/queries.txt > $BIN/b.txt 2> /dev/null
$PROG query -F -k 3 -p $TESTS/integration_b.fa -q $TESTS/queries.txt -f xor > $BIN/b_xor.txt 2> /dev/null
$PROG query -F -k 3 -p $BIN/merged.fa -q $TESTS/queries.txt > $BIN/merged.txt 2> /dev/null
$PROG query -F -k 3 -p $BIN/merged.fa -q $TESTS/queries.txt -f xor > $BIN/merged_xor.txt 2> /dev/null

#$PROG normalize -k 3 -p $BIN/merged.fa -s -l > $BIN/merged_normalized.fa 2> /dev/null
$PROG normalize -k 3 -p $BIN/merged.fa -s > $BIN/merged_normalized2.fa 2> /dev/null



diff $TESTS/result_a_complements.txt $BIN/a.txt || exit 1
echo "a.txt OK"
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

