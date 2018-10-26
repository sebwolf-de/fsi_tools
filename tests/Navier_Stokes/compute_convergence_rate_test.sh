#!/bin/sh

if [[ $1 == *"1"* ]]; then
  echo "Calculate BDF1 solutions."
fi
if [[ $1 == *"2"* ]]; then
  echo "Calculate BDF2 solutions."
fi
if [[ $1 == *"T"* ]]; then
  echo "Calculate Theta solutions."
fi
if [[ $1 == *"R"* ]]; then
  echo "Calculate reference solutions."
fi


#compute BDF1 solutions

if [[ $1 == *"1"* ]]; then
  python2 stokes_test.py convergence/test/BDF1_1.json
  cp results/ns_test_32_BDF1_dt25em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Test/BDF1_1

  python2 stokes_test.py convergence/test/BDF1_2.json
  cp results/ns_test_32_BDF1_dt125em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Test/BDF1_2

  python2 stokes_test.py convergence/test/BDF1_3.json
  cp results/ns_test_32_BDF1_dt625em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Test/BDF1_3

  python2 stokes_test.py convergence/test/BDF1_4.json
  cp results/ns_test_32_BDF1_dt3125em4_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Test/BDF1_4

  python2 stokes_test.py convergence/test/BDF1_5.json
  cp results/ns_test_32_BDF1_d15625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Test/BDF1_5

  python2 stokes_test.py convergence/test/BDF1_6.json
  cp results/ns_test_32_BDF1_dt78125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Test/BDF1_6
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 stokes_test.py convergence/test/BDF2_1
  cp results/ns_test_32_BDF2_dt1em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Test/BDF2_1

  python2 stokes_test.py convergence/test/BDF2_2
  cp results/ns_test_32_BDF2_dt125em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Test/BDF2_2

  python2 stokes_test.py convergence/test/BDF2_3
  cp results/ns_test_32_BDF2_dt625em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Test/BDF2_3

  python2 stokes_test.py convergence/test/BDF2_4.json
  cp results/ns_test_32_BDF2_dt3125em4_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Test/BDF2_4

  python2 stokes_test.py convergence/test/BDF2_5.json
  cp results/ns_test_32_BDF2_dt15625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Test/BDF2_5

  python2 stokes_test.py convergence/test/BDF2_6.json
  cp results/ns_test_32_BDF2_dt78125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Test/BDF2_6
fi
#
# #compute Theta solutions
#
# if [[ $1 == *"T"* ]]; then
#   python2 stokes_test.py convergence/test/Theta_1
#   cp results/ns_test_32_Theta_dt1em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Test/Theta_dt=1_2_result
#
#   python2 stokes_test.py convergence/test/Theta_2
#   cp results/ns_test_32_Theta_dt125em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Test/Theta_dt=1_4_result
#
#   python2 stokes_test.py convergence/test/Theta_3
#   cp results/ns_test_32_Theta_dt625em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Test/Theta_dt=1_8_result
#
#   python2 stokes_test.py convergence/test/Theta_4.json
#   cp results/ns_test_32_Theta_dt3125em4_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Test/Theta_dt=1_16_result
#
#   python2 stokes_test.py convergence/test/Theta_5.json
#   cp results/ns_test_32_Theta_dt15625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Test/Theta_dt=1_32_result
#
#   python2 stokes_test.py convergence/test/Theta_6.json
#   cp results/ns_test_32_Theta_dt78125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Test/Theta_dt=1_64_result
# fi

#compute reference solutions
cp results/ns_test_32_BDF2_dt78125em7_hx32_nu1/binary_data/mesh results/Convergence_Analysis_Test/mesh
cp results/ns_test_32_BDF2_dt78125em7_hx32_nu1/binary_data/analytical results/Convergence_Analysis_Test/reference

#compare the convergence/test rates

python2 compare.py Test
