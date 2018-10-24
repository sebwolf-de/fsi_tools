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
  python2 dlm_master.py convergence/channel/BDF1_1.json
  cp results/dlm_channel_16_BDF1_dt1em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Channel/BDF1_1

  python2 dlm_master.py convergence/channel/BDF1_2.json
  cp results/dlm_channel_16_BDF1_dt5em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/BDF1_2

  python2 dlm_master.py convergence/channel/BDF1_3.json
  cp results/dlm_channel_16_BDF1_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/BDF1_3

  python2 dlm_master.py convergence/channel/BDF1_4
  cp results/dlm_channel_16_BDF1_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/BDF1_4

  python2 dlm_master.py convergence/channel/BDF1_5.json
  cp results/dlm_channel_16_BDF1_dt625em5_hx16_hs16_k50_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Channel/BDF1_5

  python2 dlm_master.py convergence/channel/BDF1_6.json
  cp results/dlm_channel_16_BDF1_dt3125em6_hx16_hs16_k50_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Channel/BDF1_6
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 dlm_master.py convergence/channel/BDF2_1.json
  cp results/dlm_channel_16_BDF2_dt1em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Channel/BDF2_1

  python2 dlm_master.py convergence/channel/BDF2_2.json
  cp results/dlm_channel_16_BDF2_dt5em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/BDF2_2

  python2 dlm_master.py convergence/channel/BDF2_3.json
  cp results/dlm_channel_16_BDF2_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/BDF2_3

  python2 dlm_master.py convergence/channel/BDF2_4.json
  cp results/dlm_channel_16_BDF2_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/BDF2_4

  python2 dlm_master.py convergence/channel/BDF2_5.json
  cp results/dlm_channel_16_BDF2_dt625em5_hx16_hs16_k50_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Channel/BDF2_5

  python2 dlm_master.py convergence/channel/BDF2_6.json
  cp results/dlm_channel_16_BDF2_dt3125em6_hx16_hs16_k50_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Channel/BDF2_6
fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 dlm_master.py convergence/channel/Theta_1.json
  cp results/dlm_channel_16_Theta_dt1em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Channel/Theta_1

  python2 dlm_master.py convergence/channel/Theta_2.json
  cp results/dlm_channel_16_Theta_dt5em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/Theta_2

  python2 dlm_master.py convergence/channel/Theta_3.json
  cp results/dlm_channel_16_Theta_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/Theta_3

  python2 dlm_master.py convergence/channel/Theta_4.json
  cp results/dlm_channel_16_Theta_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/Theta_4

  python2 dlm_master.py convergence/channel/BDF1_5.json
  cp results/dlm_channel_16_Theta_dt625em5_hx16_hs16_k50_nu0.1/binary_data/cn_time_016 results/Convergence_Analysis_Channel/Theta_5

  python2 dlm_master.py convergence/channel/Theta_6.json
  cp results/dlm_channel_16_Theta_dt3125em6_hx16_hs16_k50_nu0.1/binary_data/cn_time_032 results/Convergence_Analysis_Channel/Theta_6
fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 dlm_master.py convergence/channel/ref.json
  cp results/dlm_channel_16_BDF1_dt5em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_200 results/Convergence_Analysis_Channel/reference
  cp results/dlm_channel_16_BDF1_dt5em4_hx16_hs16_k50_nu0.1/binary_data/mesh results/Convergence_Analysis_Channel/mesh
fi

#compare the convergence rates

python2 compare.py Channel
