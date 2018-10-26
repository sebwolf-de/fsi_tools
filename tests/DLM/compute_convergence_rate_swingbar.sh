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
  python2 dlm_master.py convergence/swingbar/BDF1_1.json
  cp results/dlm_swingbar_16_BDF1_dt2em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Swingbar/BDF1_1

  python2 dlm_master.py convergence/swingbar/BDF1_2.json
  cp results/dlm_swingbar_16_BDF1_dt1em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Swingbar/BDF1_2

  python2 dlm_master.py convergence/swingbar/BDF1_3.json
  cp results/dlm_swingbar_16_BDF1_dt5em1_hx16_hs16_k10_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Swingbar/BDF1_3

  python2 dlm_master.py convergence/swingbar/BDF1_4.json
  cp results/dlm_swingbar_16_BDF1_dt25em2_hx16_hs16_k10_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Swingbar/BDF1_4

  python2 dlm_master.py convergence/swingbar/BDF1_5.json
  cp results/dlm_swingbar_16_BDF1_dt125em3_hx16_hs16_k10_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Swingbar/BDF1_5

  python2 dlm_master.py convergence/swingbar/BDF1_6.json
  cp results/dlm_swingbar_16_BDF1_dt625em4_hx16_hs16_k10_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Swingbar/BDF1_6
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 dlm_master.py convergence/swingbar/BDF2_1.json
  cp results/dlm_swingbar_16_BDF2_dt2em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Swingbar/BDF2_1

  python2 dlm_master.py convergence/swingbar/BDF2_2.json
  cp results/dlm_swingbar_16_BDF2_dt1em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Swingbar/BDF2_2

  python2 dlm_master.py convergence/swingbar/BDF2_3.json
  cp results/dlm_swingbar_16_BDF2_dt5em1_hx16_hs16_k10_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Swingbar/BDF2_3

  python2 dlm_master.py convergence/swingbar/BDF2_4.json
  cp results/dlm_swingbar_16_BDF2_dt25em2_hx16_hs16_k10_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Swingbar/BDF2_4

  python2 dlm_master.py convergence/swingbar/BDF2_5.json
  cp results/dlm_swingbar_16_BDF2_dt125em3_hx16_hs16_k10_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Swingbar/BDF2_5

  python2 dlm_master.py convergence/swingbar/BDF2_6.json
  cp results/dlm_swingbar_16_BDF2_dt625em4_hx16_hs16_k10_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Swingbar/BDF2_6
fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 dlm_master.py convergence/swingbar/Theta_1.json
  cp results/dlm_swingbar_16_Theta_dt2em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Swingbar/Theta_1

  python2 dlm_master.py convergence/swingbar/Theta_2.json
  cp results/dlm_swingbar_16_Theta_dt1em0_hx16_hs16_k10_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Swingbar/Theta_2

  python2 dlm_master.py convergence/swingbar/Theta_3.json
  cp results/dlm_swingbar_16_Theta_dt5em1_hx16_hs16_k10_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Swingbar/Theta_3

  python2 dlm_master.py convergence/swingbar/Theta_4.json
  cp results/dlm_swingbar_16_Theta_dt25em2_hx16_hs16_k10_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Swingbar/Theta_4

  python2 dlm_master.py convergence/swingbar/Theta_5.json
  cp results/dlm_swingbar_16_Theta_dt125em3_hx16_hs16_k10_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Swingbar/Theta_5

  python2 dlm_master.py convergence/swingbar/Theta_6.json
  cp results/dlm_swingbar_16_Theta_dt625em4_hx16_hs16_k10_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Swingbar/Theta_6
fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 dlm_master.py convergence/swingbar/ref.json
  cp results/dlm_swingbar_16_BDF1_dt5em3_hx16_hs16_k10_nu0.1/binary_data/cn_time_200 results/Convergence_Analysis_Swingbar/reference
  cp results/dlm_swingbar_16_BDF1_dt5em3_hx16_hs16_k10_nu0.1/binary_data/mesh results/Convergence_Analysis_Swingbar/mesh
fi

#compare the convergence/swingbar rates

python2 compare.py Swingbar
