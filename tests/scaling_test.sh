#!/bin/bash -x

echo > scaling_test.output

./spamm_multiply_01 >> scaling_test.output
./spamm_multiply_02 >> scaling_test.output
./spamm_multiply_04 >> scaling_test.output
./spamm_multiply_08 >> scaling_test.output
./spamm_multiply_12 >> scaling_test.output
./spamm_multiply_16 >> scaling_test.output
./spamm_multiply_20 >> scaling_test.output
./spamm_multiply_24 >> scaling_test.output
./spamm_multiply_28 >> scaling_test.output
./spamm_multiply_32 >> scaling_test.output
./spamm_multiply_36 >> scaling_test.output
./spamm_multiply_40 >> scaling_test.output
./spamm_multiply_44 >> scaling_test.output
./spamm_multiply_48 >> scaling_test.output

grep SpAMM_SCALING scaling_test.output | awk '{print $2, $3}' > scaling_test.dat
