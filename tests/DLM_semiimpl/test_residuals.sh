#!/bin/bash

echo 'BDF1, dt=25em3'
python3 dlm_semiimplicit.py params/res/1.json | grep residual
echo 'BDF2, dt=25em3'
python3 dlm_semiimplicit.py params/res/2.json | grep residual
echo 'Theta, dt=25em3'
python3 dlm_semiimplicit.py params/res/3.json | grep residual
echo 'BDF1, dt=25em4'
python3 dlm_semiimplicit.py params/res/4.json | grep residual
echo 'BDF2, dt=25em4'
python3 dlm_semiimplicit.py params/res/5.json | grep residual
echo 'Theta, dt=25em4'
python3 dlm_semiimplicit.py params/res/6.json | grep residual
