#!/bin/sh

if [[ $1 == *"1"* ]]; then
  echo "Calculate BDF1 solutions."
fi
if [[ $1 == *"2"* ]]; then
  echo "Calculate BDF2 solutions."
fi
if [[ $1 == *"R"* ]]; then
  echo "Calculate reference solutions."
fi


#compute BDF1 solutions

if [[ $1 == *"1"* ]]; then
  python2 ns_master.py convergence/BDF1_dt\=1_2.json
  cp results/ns_cavity_32_BDF1_dt1em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Cavity/BDF1_dt=1_2_result

  python2 ns_master.py convergence/BDF1_dt\=1_4.json
  cp results/ns_cavity_32_BDF1_dt5em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF1_dt=1_4_result

  python2 ns_master.py convergence/BDF1_dt\=1_8.json
  cp results/ns_cavity_32_BDF1_dt25em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF1_dt=1_8_result

  python2 ns_master.py convergence/BDF1_dt\=1_16.json
  cp results/ns_cavity_32_BDF1_dt125em5_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF1_dt=1_16_result

  python2 ns_master.py convergence/BDF1_dt\=1_32.json
  cp results/ns_cavity_32_BDF1_d625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF1_dt=1_32_result

  python2 ns_master.py convergence/BDF1_dt\=1_64.json
  cp results/ns_cavity_32_BDF1_dt3125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF1_dt=1_64_result
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 ns_master.py convergence/BDF2_dt\=1_2.json
  cp results/ns_cavity_32_BDF2_dt1em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Cavity/BDF2_dt=1_2_result

  python2 ns_master.py convergence/BDF2_dt\=1_4.json
  cp results/ns_cavity_32_BDF2_dt5em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF2_dt=1_4_result

  python2 ns_master.py convergence/BDF2_dt\=1_8.json
  cp results/ns_cavity_32_BDF2_dt25em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF2_dt=1_8_result

  python2 ns_master.py convergence/BDF2_dt\=1_16.json
  cp results/ns_cavity_32_BDF2_dt125em5_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF2_dt=1_16_result

  python2 ns_master.py convergence/BDF2_dt\=1_32.json
  cp results/ns_cavity_32_BDF2_dt625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF2_dt=1_32_result

  python2 ns_master.py convergence/BDF2_dt\=1_64.json
  cp results/ns_cavity_32_BDF2_dt3125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF2_dt=1_64_result
fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 ns_master.py convergence/Theta_dt\=1_2.json
  cp results/ns_cavity_32_Theta_dt1em2_hx32_nu1/binary_data/cn_time_001 results/Convergence_Analysis_Cavity/Theta_dt=1_2_result

  python2 ns_master.py convergence/Theta_dt\=1_4.json
  cp results/ns_cavity_32_Theta_dt5em3_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/Theta_dt=1_4_result

  python2 ns_master.py convergence/Theta_dt\=1_8.json
  cp results/ns_cavity_32_Theta_dt25em4_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/Theta_dt=1_8_result

  python2 ns_master.py convergence/Theta_dt\=1_16.json
  cp results/ns_cavity_32_Theta_dt125em5_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/Theta_dt=1_16_result

  python2 ns_master.py convergence/Theta_dt\=1_32.json
  cp results/ns_cavity_32_Theta_dt625em6_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/Theta_dt=1_32_result

  python2 ns_master.py convergence/Theta_dt\=1_64.json
  cp results/ns_cavity_32_Theta_dt3125em7_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/Theta_dt=1_64_result
fi
#compute reference solutions

if [[ $1 == *"R"* ]]; then
python2 ns_master.py convergence/ref.json
cp results/ns_cavity_32_BDF1_dt5em5_hx32_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Cavity/reference
cp results/ns_cavity_32_BDF1_dt5em5_hx32_nu1/binary_data/mesh results/Convergence_Analysis_Cavity/mesh
fi

python2 compare.py

#compare the convergence rates

python2 compare.py
