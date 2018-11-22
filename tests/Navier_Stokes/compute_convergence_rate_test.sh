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
  python2 ns_master.py convergence/test/BDF1_1.json
  cp results/ns_test_12_BDF1_dt1em3_hx12_nu1/binary_data/cn_time_0002 results/Convergence_Analysis_Test/BDF1_1

  python2 ns_master.py convergence/test/BDF1_2.json
  cp results/ns_test_12_BDF1_dt5em4_hx12_nu1/binary_data/cn_time_0004 results/Convergence_Analysis_Test/BDF1_2

  python2 ns_master.py convergence/test/BDF1_3.json
  cp results/ns_test_12_BDF1_dt25em5_hx12_nu1/binary_data/cn_time_0008 results/Convergence_Analysis_Test/BDF1_3

  python2 ns_master.py convergence/test/BDF1_4.json
  cp results/ns_test_12_BDF1_dt125em6_hx12_nu1/binary_data/cn_time_0016 results/Convergence_Analysis_Test/BDF1_4

fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 ns_master.py convergence/test/BDF2_1.json
  cp results/ns_test_12_BDF2_dt1em3_hx12_nu1/binary_data/cn_time_0002 results/Convergence_Analysis_Test/BDF2_1

  python2 ns_master.py convergence/test/BDF2_2.json
  cp results/ns_test_12_BDF2_dt5em4_hx12_nu1/binary_data/cn_time_0004 results/Convergence_Analysis_Test/BDF2_2

  python2 ns_master.py convergence/test/BDF2_3.json
  cp results/ns_test_12_BDF2_dt25em5_hx12_nu1/binary_data/cn_time_0008 results/Convergence_Analysis_Test/BDF2_3

  python2 ns_master.py convergence/test/BDF2_4.json
  cp results/ns_test_12_BDF2_dt125em6_hx12_nu1/binary_data/cn_time_0016 results/Convergence_Analysis_Test/BDF2_4
fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
    python2 ns_master.py convergence/test/Theta_1.json
    cp results/ns_test_12_Theta_dt1em3_hx12_nu1/binary_data/cn_time_0002 results/Convergence_Analysis_Test/Theta_1

    python2 ns_master.py convergence/test/Theta_2.json
    cp results/ns_test_12_Theta_dt5em4_hx12_nu1/binary_data/cn_time_0004 results/Convergence_Analysis_Test/Theta_2

    python2 ns_master.py convergence/test/Theta_3.json
    cp results/ns_test_12_Theta_dt25em5_hx12_nu1/binary_data/cn_time_0008 results/Convergence_Analysis_Test/Theta_3

    python2 ns_master.py convergence/test/Theta_4.json
    cp results/ns_test_12_Theta_dt125em6_hx12_nu1/binary_data/cn_time_0016 results/Convergence_Analysis_Test/Theta_4
fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 ns_master.py convergence/test/ref.json
  cp results/ns_test_12_BDF2_dt1em5_hx12_nu1/binary_data/mesh results/Convergence_Analysis_Test/mesh
  cp results/ns_test_12_BDF2_dt1em5_hx12_nu1/binary_data/cn_time_0200 results/Convergence_Analysis_Test/reference
fi

#compare the convergence/test rates

python2 compare.py Test
