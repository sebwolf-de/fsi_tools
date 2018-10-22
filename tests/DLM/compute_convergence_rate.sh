#!/bin/sh

#compute BDF1 solutions
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

if [[ $1 == *"1"* ]]; then
  python2 dlm_master.py convergence/BDF1_dt\=1_2.json
  cp results/dlm_annulus_8_BDF1_dt5em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF1_dt=1_2_result
  cp results/dlm_annulus_8_BDF1_dt5em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_2_time

  python2 dlm_master.py convergence/BDF1_dt\=1_4.json
  cp results/dlm_annulus_8_BDF1_dt25em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF1_dt=1_4_result
  cp results/dlm_annulus_8_BDF1_dt25em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_4_time

  python2 dlm_master.py convergence/BDF1_dt\=1_8.json
  cp results/dlm_annulus_8_BDF1_dt125em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF1_dt=1_8_result
  cp results/dlm_annulus_8_BDF1_dt125em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_8_time

  python2 dlm_master.py convergence/BDF1_dt\=1_16.json
  cp results/dlm_annulus_8_BDF1_dt625em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF1_dt=1_16_result
  cp results/dlm_annulus_8_BDF1_dt625em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_16_time

  python2 dlm_master.py convergence/BDF1_dt\=1_32.json
  cp results/dlm_annulus_8_BDF1_dt3125em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF1_dt=1_32_result
  cp results/dlm_annulus_8_BDF1_dt3125em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_32_time

  python2 dlm_master.py convergence/BDF1_dt\=1_64.json
  cp results/dlm_annulus_8_BDF1_dt15625em6_hx8_hs8_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/BDF1_dt=1_64_result
  cp results/dlm_annulus_8_BDF1_dt15625em6_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_dt=1_64_time
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 dlm_master.py convergence/BDF2_dt\=1_2.json
  cp results/dlm_annulus_8_BDF2_dt5em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF2_dt=1_2_result
  cp results/dlm_annulus_8_BDF2_dt5em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_2_time

  python2 dlm_master.py convergence/BDF2_dt\=1_4.json
  cp results/dlm_annulus_8_BDF2_dt25em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF2_dt=1_4_result
  cp results/dlm_annulus_8_BDF2_dt25em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_4_time

  python2 dlm_master.py convergence/BDF2_dt\=1_8.json
  cp results/dlm_annulus_8_BDF2_dt125em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF2_dt=1_8_result
  cp results/dlm_annulus_8_BDF2_dt125em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_8_time

  python2 dlm_master.py convergence/BDF2_dt\=1_16.json
  cp results/dlm_annulus_8_BDF2_dt625em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF2_dt=1_16_result
  cp results/dlm_annulus_8_BDF2_dt625em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_16_time

  python2 dlm_master.py convergence/BDF2_dt\=1_32.json
  cp results/dlm_annulus_8_BDF2_dt3125em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF2_dt=1_32_result
  cp results/dlm_annulus_8_BDF2_dt3125em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_32_time

  python2 dlm_master.py convergence/BDF2_dt\=1_64.json
  cp results/dlm_annulus_8_BDF2_dt15625em6_hx8_hs8_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/BDF2_dt=1_64_result
  cp results/dlm_annulus_8_BDF2_dt15625em6_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_dt=1_64_time
fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 dlm_master.py convergence/Theta_dt\=1_2.json
  cp results/dlm_annulus_8_Theta_dt5em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/Theta_dt=1_2_result
  cp results/dlm_annulus_8_Theta_dt5em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_2_time

  python2 dlm_master.py convergence/Theta_dt\=1_4.json
  cp results/dlm_annulus_8_Theta_dt25em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/Theta_dt=1_4_result
  cp results/dlm_annulus_8_Theta_dt25em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_4_time

  python2 dlm_master.py convergence/Theta_dt\=1_8.json
  cp results/dlm_annulus_8_Theta_dt125em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/Theta_dt=1_8_result
  cp results/dlm_annulus_8_Theta_dt125em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_8_time

  python2 dlm_master.py convergence/Theta_dt\=1_16.json
  cp results/dlm_annulus_8_Theta_dt625em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/Theta_dt=1_16_result
  cp results/dlm_annulus_8_Theta_dt625em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_16_time

  python2 dlm_master.py convergence/Theta_dt\=1_32.json
  cp results/dlm_annulus_8_Theta_dt3125em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/Theta_dt=1_32_result
  cp results/dlm_annulus_8_Theta_dt3125em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_32_time

  python2 dlm_master.py convergence/Theta_dt\=1_64.json
  cp results/dlm_annulus_8_Theta_dt15625em6_hx8_hs8_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/Theta_dt=1_64_result
  cp results/dlm_annulus_8_Theta_dt15625em6_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/Theta_dt=1_64_time
fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 dlm_master.py convergence/ref.json
  cp results/dlm_annulus_8_BDF1_dt5em3_hx8_hs8_k10_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Annulus/reference
  cp results/dlm_annulus_8_BDF1_dt5em3_hx8_hs8_k10_nu1/binary_data/mesh results/Convergence_Analysis_Annulus/mesh
fi

#compare the convergence rates

python2 compare.py
