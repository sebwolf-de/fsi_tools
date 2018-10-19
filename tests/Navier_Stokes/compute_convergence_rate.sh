#!/bin/sh

#compute BDF1 solutions

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_2.json
#cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF1_dt=1_2_result
#cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_2_time

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_4.json
#cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF1_dt=1_4_result
#cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_4_time

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_8.json
#cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF1_dt=1_8_result
#cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_8_time

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_16.json
#cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF1_dt=1_16_result
#cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_16_time

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_32.json
#cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF1_dt=1_32_result
#cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_32_time

#python2 ns_lid_cavity.py convergence/BDF1_dt\=1_64.json
#cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/BDF1_dt=1_64_result
#cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF1_dt=1_64_time

#compute BDF2 solutions

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_2.json
cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/BDF2_dt=1_2_result
#cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_2_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_4.json
cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/BDF2_dt=1_4_result
#cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_4_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_8.json
cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/BDF2_dt=1_8_result
#cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_8_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_16.json
cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/BDF2_dt=1_16_result
#cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_16_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_32.json
cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/BDF2_dt=1_32_result
#cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_32_time

python2 ns_lid_cavity.py convergence/BDF2_dt\=1_64.json
cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/BDF2_dt=1_64_result
#cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/BDF2_dt=1_64_time

#compute Theta solutions

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_2.json
#cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_Cavity/Theta_dt=1_2_result
#cp results/ns_cavity_32_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_2_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_4.json
#cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_Cavity/Theta_dt=1_4_result
#cp results/ns_cavity_32_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_4_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_8.json
#cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_Cavity/Theta_dt=1_8_result
#cp results/ns_cavity_32_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_8_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_16.json
#cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_Cavity/Theta_dt=1_16_result
#cp results/ns_cavity_32_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_16_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_32.json
#cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_Cavity/Theta_dt=1_32_result
#cp results/ns_cavity_32_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_32_time

#python2 ns_lid_cavity.py convergence/Theta_dt\=1_64.json
#cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_Cavity/Theta_dt=1_64_result
#cp results/ns_cavity_32_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_Cavity/Theta_dt=1_64_time

#compute reference solutions

python2 ns_lid_cavity.py convergence/ref.json
cp results/ns_cavity_32_dt5em3_hx32_nu1/binary_data/cn_time_200 results/Convergence_Analysis_Cavity/reference
#cp results/ns_cavity_32_dt5em3_hx32_nu1/binary_data/mesh results/Convergence_Analysis_Cavity/mesh

#compare the convergence rates

python2 compare.py
