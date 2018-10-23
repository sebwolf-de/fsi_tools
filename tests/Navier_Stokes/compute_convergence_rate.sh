#!/bin/sh

#compute BDF1 solutions

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_2.json
cp results/ns_cavity_8_BDF1_dt5em1_hx8_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF1_dt=1_2_result
#cp results/ns_cavity_8_dt5em1_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_2_time

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_4.json
cp results/ns_cavity_8_BDF1_dt25em2_hx8_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF1_dt=1_4_result
#cp results/ns_cavity_8_dt25em2_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_4_time

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_8.json
cp results/ns_cavity_8_BDF1_dt125em3_hx8_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF1_dt=1_8_result
#cp results/ns_cavity_8_dt125em3_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_8_time

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_16.json
cp results/ns_cavity_8_BDF1_dt625em4_hx8_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF1_dt=1_16_result
#cp results/ns_cavity_8_dt625em4_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_16_time

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_32.json
cp results/ns_cavity_8_BDF1_dt3125em5_hx8_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF1_dt=1_32_result
#cp results/ns_cavity_8_dt3125em5_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_32_time

python2 ns_lid_cavity.py convergence/BDF1_dt\=1_64.json
cp results/ns_cavity_8_BDF1_dt15625em6_hx8_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/BDF1_dt=1_64_result
#cp results/ns_cavity_8_dt15625em6_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_64_time

#compute BDF2 solutions

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_2.json
cp results/ns_cavity_8_BDF2_dt5em1_hx8_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF2_dt=1_2_result
#cp results/ns_cavity_8_dt5em1_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_2_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_4.json
cp results/ns_cavity_8_BDF2_dt25em2_hx8_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF2_dt=1_4_result
#cp results/ns_cavity_8_dt25em2_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_4_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_8.json
cp results/ns_cavity_8_BDF2_dt125em3_hx8_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF2_dt=1_8_result
#cp results/ns_cavity_8_dt125em3_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_8_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_16.json
cp results/ns_cavity_8_BDF2_dt625em4_hx8_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF2_dt=1_16_result
#cp results/ns_cavity_8_dt625em4_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_16_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_32.json
cp results/ns_cavity_8_BDF2_dt3125em5_hx8_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF2_dt=1_32_result
#cp results/ns_cavity_8_dt3125em5_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_32_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_64.json
cp results/ns_cavity_8_BDF2_dt15625em6_hx8_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/BDF2_dt=1_64_result
#cp results/ns_cavity_8_dt15625em6_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_64_time

#compute Theta solutions

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_2.json
#cp results/ns_cavity_8_Theta_dt5em1_hx8_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/Theta_dt=1_2_result
#cp results/ns_cavity_8_dt5em1_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_2_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_4.json
#cp results/ns_cavity_8_Theta_dt25em2_hx8_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/Theta_dt=1_4_result
#cp results/ns_cavity_8_dt25em2_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_4_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_8.json
#cp results/ns_cavity_8_Theta_dt125em3_hx8_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/Theta_dt=1_8_result
#cp results/ns_cavity_8_dt125em3_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_8_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_16.json
#cp results/ns_cavity_8_Theta_dt625em4_hx8_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/Theta_dt=1_16_result
#cp results/ns_cavity_8_dt625em4_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_16_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_32.json
#cp results/ns_cavity_8_Theta_dt3125em5_hx8_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/Theta_dt=1_32_result
#cp results/ns_cavity_8_dt3125em5_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_32_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_64.json
#cp results/ns_cavity_8_Theta_dt15625em6_hx8_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/Theta_dt=1_64_result
#cp results/ns_cavity_8_dt15625em6_hx8_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_64_time

#compute reference solutions

python2 ns_lid_cavity.py convergence/ref.json
cp results/ns_cavity_8_BDF1_dt5em3_hx8_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Cavity/reference
cp results/ns_cavity_8_BDF1_dt5em3_hx8_nu1/binary_data/mesh results/Convergence_Analysis_Cavity/mesh

#compare the convergence rates

python2 compare.py
