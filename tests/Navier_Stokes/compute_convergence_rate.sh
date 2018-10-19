#!/bin/sh

#compute BDF1 solutions

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_2.json
#cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_NS/BDF1_dt=1_2_result
#cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_2_time

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_4.json
#cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_NS/BDF1_dt=1_4_result
#cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_4_time

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_8.json
#cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_NS/BDF1_dt=1_8_result
#cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_8_time

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_16.json
#cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_NS/BDF1_dt=1_16_result
#cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_16_time

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_32.json
#cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_NS/BDF1_dt=1_32_result
#cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_32_time

#python ns_lid_cavity.py convergence_ns/BDF1_dt\=1_64.json
#cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_NS/BDF1_dt=1_64_result
#cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF1_dt=1_64_time

#compute BDF2 solutions

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_2.json
cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_NS/BDF2_dt=1_2_result
cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_2_time

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_4.json
cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_NS/BDF2_dt=1_4_result
cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_4_time

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_8.json
cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_NS/BDF2_dt=1_8_result
cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_8_time

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_16.json
cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_NS/BDF2_dt=1_16_result
cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_16_time

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_32.json
cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_NS/BDF2_dt=1_32_result
cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_32_time

python ns_lid_cavity.py convergence_ns/BDF2_dt\=1_64.json
cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_NS/BDF2_dt=1_64_result
cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/BDF2_dt=1_64_time

#compute Theta solutions

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_2.json
#cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/cn_time_002 results/Convergence_Analysis_NS/Theta_dt=1_2_result
#cp results/ns_cavity_8_dt5em1_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_2_time

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_4.json
#cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/cn_time_004 results/Convergence_Analysis_NS/Theta_dt=1_4_result
#cp results/ns_cavity_8_dt25em2_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_4_time

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_8.json
#cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/cn_time_008 results/Convergence_Analysis_NS/Theta_dt=1_8_result
#cp results/ns_cavity_8_dt125em3_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_8_time

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_16.json
#cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/cn_time_016 results/Convergence_Analysis_NS/Theta_dt=1_16_result
#cp results/ns_cavity_8_dt625em4_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_16_time

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_32.json
#cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/cn_time_032 results/Convergence_Analysis_NS/Theta_dt=1_32_result
#cp results/ns_cavity_8_dt3125em5_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_32_time

#python ns_lid_cavity.py convergence_ns/Theta_dt\=1_64.json
#cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/cn_time_064 results/Convergence_Analysis_NS/Theta_dt=1_64_result
#cp results/ns_cavity_8_dt15625em6_hx32_nu1/binary_data/time results/Convergence_Analysis_NS/Theta_dt=1_64_time

#compute reference solutions

python ns_lid_cavity.py convergence_ns/ref.json
cp results/ns_cavity_8_dt5em3_hx32_nu1/binary_data/cn_time_200 results/Convergence_Analysis_NS/reference
cp results/ns_cavity_8_dt5em3_hx32_nu1/binary_data/mesh results/Convergence_Analysis_NS/mesh

#compare the convergence_ns rates

python compare.py
