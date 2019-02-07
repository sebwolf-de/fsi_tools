#!/bin/sh

if [[ $1 == *"1"* ]]; then
  echo "Calculate BDF1 solutions."
fi
if [[ $1 == *"2"* ]]; then
  echo "Calculate BDF2 solutions."
fi
if [[ $1 == *"C"* ]]; then
  echo "Calculate Crank-Nicolson solutions."
fi
if [[ $1 == *"T"* ]]; then
  echo "Calculate Trapezoidal solutions."
fi
if [[ $1 == *"R"* ]]; then
  echo "Calculate reference solutions."
fi

mkdir -p results/Convergence_Analysis_Annulus

#compute BDF1 solutions

if [[ $1 == *"1"* ]]; then
  python3 dlm_semiimplicit.py params/conv/BDF1_dt\=1_2.json
  python3 dlm_semiimplicit.py params/conv/BDF1_dt\=1_4.json
  python3 dlm_semiimplicit.py params/conv/BDF1_dt\=1_8.json
  python3 dlm_semiimplicit.py params/conv/BDF1_dt\=1_16.json
  python3 dlm_semiimplicit.py params/conv/BDF1_dt\=1_32.json
fi
cp results/dlm_annulus_8_BDF1_dt1em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF1_1
cp results/dlm_annulus_8_BDF1_dt1em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_1_time

cp results/dlm_annulus_8_BDF1_dt5em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF1_2
cp results/dlm_annulus_8_BDF1_dt5em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_2_time

cp results/dlm_annulus_8_BDF1_dt25em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF1_3
cp results/dlm_annulus_8_BDF1_dt25em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_3_time

cp results/dlm_annulus_8_BDF1_dt125em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF1_4
cp results/dlm_annulus_8_BDF1_dt125em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_4_time

cp results/dlm_annulus_8_BDF1_dt625em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF1_5
cp results/dlm_annulus_8_BDF1_dt625em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF1_5_time

#compute BDF2 solutions

if [[ $1 == *"2"* ]]; then
  python3 dlm_semiimplicit.py params/conv/BDF2_dt\=1_2.json
  python3 dlm_semiimplicit.py params/conv/BDF2_dt\=1_4.json
  python3 dlm_semiimplicit.py params/conv/BDF2_dt\=1_8.json
  python3 dlm_semiimplicit.py params/conv/BDF2_dt\=1_16.json
  python3 dlm_semiimplicit.py params/conv/BDF2_dt\=1_32.json
fi

cp results/dlm_annulus_8_BDF2_dt1em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/BDF2_1
cp results/dlm_annulus_8_BDF2_dt1em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_1_time

cp results/dlm_annulus_8_BDF2_dt5em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/BDF2_2
cp results/dlm_annulus_8_BDF2_dt5em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_2_time

cp results/dlm_annulus_8_BDF2_dt25em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/BDF2_3
cp results/dlm_annulus_8_BDF2_dt25em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_3_time

cp results/dlm_annulus_8_BDF2_dt125em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/BDF2_4
cp results/dlm_annulus_8_BDF2_dt125em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_4_time

cp results/dlm_annulus_8_BDF2_dt625em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/BDF2_5
cp results/dlm_annulus_8_BDF2_dt625em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/BDF2_5_time

#compute CN solutions

if [[ $1 == *"C"* ]]; then
  python3 dlm_semiimplicit.py params/conv/CN_dt\=1_2.json
  python3 dlm_semiimplicit.py params/conv/CN_dt\=1_4.json
  python3 dlm_semiimplicit.py params/conv/CN_dt\=1_8.json
  python3 dlm_semiimplicit.py params/conv/CN_dt\=1_16.json
  python3 dlm_semiimplicit.py params/conv/CN_dt\=1_32.json
fi

cp results/dlm_annulus_8_CN_dt1em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/CN_1
cp results/dlm_annulus_8_CN_dt1em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/CN_1_time

cp results/dlm_annulus_8_CN_dt5em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/CN_2
cp results/dlm_annulus_8_CN_dt5em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/CN_2_time

cp results/dlm_annulus_8_CN_dt25em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/CN_3
cp results/dlm_annulus_8_CN_dt25em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/CN_3_time

cp results/dlm_annulus_8_CN_dt125em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/CN_4
cp results/dlm_annulus_8_CN_dt125em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/CN_4_time

cp results/dlm_annulus_8_CN_dt625em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/CN_5
cp results/dlm_annulus_8_CN_dt625em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/CN_5_time

#compute TR solutions

if [[ $1 == *"T"* ]]; then
  python3 dlm_semiimplicit.py params/conv/TR_dt\=1_2.json
  python3 dlm_semiimplicit.py params/conv/TR_dt\=1_4.json
  python3 dlm_semiimplicit.py params/conv/TR_dt\=1_8.json
  python3 dlm_semiimplicit.py params/conv/TR_dt\=1_16.json
  python3 dlm_semiimplicit.py params/conv/TR_dt\=1_32.json
fi

cp results/dlm_annulus_8_TR_dt1em1_hx8_hs8_k10_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Annulus/TR_1
cp results/dlm_annulus_8_TR_dt1em1_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/TR_1_time

cp results/dlm_annulus_8_TR_dt5em2_hx8_hs8_k10_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Annulus/TR_2
cp results/dlm_annulus_8_TR_dt5em2_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/TR_2_time

cp results/dlm_annulus_8_TR_dt25em3_hx8_hs8_k10_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Annulus/TR_3
cp results/dlm_annulus_8_TR_dt25em3_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/TR_3_time

cp results/dlm_annulus_8_TR_dt125em4_hx8_hs8_k10_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Annulus/TR_4
cp results/dlm_annulus_8_TR_dt125em4_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/TR_4_time

cp results/dlm_annulus_8_TR_dt625em5_hx8_hs8_k10_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Annulus/TR_5
cp results/dlm_annulus_8_TR_dt625em5_hx8_hs8_k10_nu1/binary_data/time results/Convergence_Analysis_Annulus/TR_5_time

#compute reference solutions

if [[ $1 == *"R"* ]]; then
  python3 dlm_semiimplicit.py params/conv/ref.json
fi

cp results/dlm_annulus_8_BDF2_dt1em3_hx8_hs8_k10_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Annulus/reference
cp results/dlm_annulus_8_BDF2_dt1em1_hx8_hs8_k10_nu1/binary_data/mesh results/Convergence_Analysis_Annulus/mesh

#compare the params/conv rates

python3 compare.py Annulus
