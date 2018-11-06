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
  python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_2.json
  cp results/dlm_annulus_16_BDF1_dt1em1_hx16_hs16_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF1_1

  python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_4.json
  cp results/dlm_annulus_16_BDF1_dt5em2_hx16_hs16_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF1_2

  python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_8.json
  cp results/dlm_annulus_16_BDF1_dt25em3_hx16_hs16_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF1_3

  python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_16.json
  cp results/dlm_annulus_16_BDF1_dt125em4_hx16_hs16_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF1_4

  #python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_32.json
  #cp results/dlm_annulus_16_BDF1_dt625em5_hx16_hs16_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF1_5

  #python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_64.json
  #cp results/dlm_annulus_16_BDF1_dt3125em6_hx16_hs16_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/BDF1_6

  #python2 dlm_fixpoint.py convergence/annulus/BDF1_dt\=1_128.json
  #cp results/dlm_annulus_16_BDF1_dt15625em7_hx16_hs16_k10_nu1/binary_data/cn_time_128 results/Convergence_Analysis_Annulus/BDF1_7
fi

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_2.json
  cp results/dlm_annulus_16_BDF2_dt1em1_hx16_hs16_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF2_1

  python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_4.json
  cp results/dlm_annulus_16_BDF2_dt5em2_hx16_hs16_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF2_2

  python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_8.json
  cp results/dlm_annulus_16_BDF2_dt25em3_hx16_hs16_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF2_3

  python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_16.json
  cp results/dlm_annulus_16_BDF2_dt125em4_hx16_hs16_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF2_4

  #python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_32.json
  #cp results/dlm_annulus_16_BDF2_dt625em5_hx16_hs16_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF2_5

  #python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_64.json
  #cp results/dlm_annulus_16_BDF2_dt3125em6_hx16_hs16_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/BDF2_6

  #python2 dlm_fixpoint.py convergence/annulus/BDF2_dt\=1_128.json
  #cp results/dlm_annulus_16_BDF2_dt15625em7_hx16_hs16_k10_nu1/binary_data/cn_time_128 results/Convergence_Analysis_Annulus/BDF2_7

fi

#compute Theta solutions

if [[ $1 == *"T"* ]]; then
  python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_2.json
  cp results/dlm_annulus_16_Theta_dt1em1_hx16_hs16_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/Theta_1

  python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_4.json
  cp results/dlm_annulus_16_Theta_dt5em2_hx16_hs16_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/Theta_2

  python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_8.json
  cp results/dlm_annulus_16_Theta_dt25em3_hx16_hs16_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/Theta_3

  python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_16.json
  cp results/dlm_annulus_16_Theta_dt125em4_hx16_hs16_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/Theta_4

  #python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_32.json
  #cp results/dlm_annulus_16_Theta_dt625em5_hx16_hs16_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/Theta_5

  #python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_64.json
  #cp results/dlm_annulus_16_Theta_dt3125em6_hx16_hs16_k10_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Annulus/Theta_6

  #python2 dlm_fixpoint.py convergence/annulus/Theta_dt\=1_128.json
  #cp results/dlm_annulus_16_Theta_dt15625em7_hx16_hs16_k10_nu1/binary_data/cn_time_128 results/Convergence_Analysis_Annulus/Theta_7
fi

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python2 dlm_fixpoint.py convergence/annulus/ref.json
  cp results/dlm_annulus_16_BDF2_dt5em4_hx16_hs16_k10_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Annulus/reference
  cp results/dlm_annulus_16_BDF2_dt5em4_hx16_hs16_k10_nu1/binary_data/mesh results/Convergence_Analysis_Annulus/mesh
fi

#compare the convergence/annulus rates

python2 compare.py Annulus
