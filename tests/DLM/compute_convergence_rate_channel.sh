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
  python2 dlm_master.py convergence/channel/BDF1_dt\=1_2.json
  cp results/dlm_channel_16_BDF1_dt1em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_001 results/Convergence_Analysis_Channel/BDF1_dt=1_2_result
  cp results/dlm_channel_16_BDF1_dt1em1_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF1_dt=1_2_time

  python2 dlm_master.py convergence/channel/BDF1_dt\=1_4.json
  cp results/dlm_channel_16_BDF1_dt5em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/BDF1_dt=1_4_result
  cp results/dlm_channel_16_BDF1_dt5em2_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF1_dt=1_4_time

  python2 dlm_master.py convergence/channel/BDF1_dt\=1_8.json
  cp results/dlm_channel_16_BDF1_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/BDF1_dt=1_8_result
  cp results/dlm_channel_16_BDF1_dt25em3_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF1_dt=1_8_time

  python2 dlm_master.py convergence/channel/BDF1_dt\=1_16.json
  cp results/dlm_channel_16_BDF1_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/BDF1_dt=1_16_result
  cp results/dlm_channel_16_BDF1_dt125em4_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF1_dt=1_16_time

fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 dlm_master.py convergence/channel/BDF2_dt\=1_2.json
  cp results/dlm_channel_16_BDF2_dt5em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/BDF2_dt=1_2_result
  cp results/dlm_channel_16_BDF2_dt5em1_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF2_dt=1_2_time

  python2 dlm_master.py convergence/channel/BDF2_dt\=1_4.json
  cp results/dlm_channel_16_BDF2_dt25em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/BDF2_dt=1_4_result
  cp results/dlm_channel_16_BDF2_dt25em2_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF2_dt=1_4_time

  python2 dlm_master.py convergence/channel/BDF2_dt\=1_8.json
  cp results/dlm_channel_16_BDF2_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/BDF2_dt=1_8_result
  cp results/dlm_channel_16_BDF2_dt25em3_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF2_dt=1_8_time

  python2 dlm_master.py convergence/channel/BDF2_dt\=1_16.json
  cp results/dlm_channel_16_BDF2_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/BDF2_dt=1_16_result
  cp results/dlm_channel_16_BDF2_dt125em4_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/BDF2_dt=1_16_time

fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 dlm_master.py convergence/channel/Theta_dt\=1_2.json
  cp results/dlm_channel_16_Theta_dt5em1_hx16_hs16_k50_nu0.1/binary_data/cn_time_002 results/Convergence_Analysis_Channel/Theta_dt=1_2_result
  cp results/dlm_channel_16_Theta_dt5em1_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/Theta_dt=1_2_time

  python2 dlm_master.py convergence/channel/Theta_dt\=1_4.json
  cp results/dlm_channel_16_Theta_dt25em2_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/Theta_dt=1_4_result
  cp results/dlm_channel_16_Theta_dt25em2_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/Theta_dt=1_4_time

  python2 dlm_master.py convergence/channel/Theta_dt\=1_8.json
  cp results/dlm_channel_16_Theta_dt25em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_004 results/Convergence_Analysis_Channel/Theta_dt=1_8_result
  cp results/dlm_channel_16_Theta_dt25em3_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/Theta_dt=1_8_time

  python2 dlm_master.py convergence/channel/Theta_dt\=1_16.json
  cp results/dlm_channel_16_Theta_dt125em4_hx16_hs16_k50_nu0.1/binary_data/cn_time_008 results/Convergence_Analysis_Channel/Theta_dt=1_16_result
  cp results/dlm_channel_16_Theta_dt125em4_hx16_hs16_k50_nu0.1/binary_data/time results/Convergence_Analysis_Channel/Theta_dt=1_16_time

fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 dlm_master.py convergence/channel/ref.json
  cp results/dlm_channel_16_BDF1_dt5em3_hx16_hs16_k50_nu0.1/binary_data/cn_time_200 results/Convergence_Analysis_Channel/reference
  cp results/dlm_channel_16_BDF1_dt5em3_hx16_hs16_k50_nu0.1/binary_data/mesh results/Convergence_Analysis_Channel/mesh
fi

#compare the convergence rates

python2 compare.py Channel
