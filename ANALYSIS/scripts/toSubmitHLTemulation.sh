##
#  script to submit jobs to perform the HLT-emulation on data
##


# 2016 PRE VFP
# 2016 POST VFP
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2016_F.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2016_G.txt

# 2017
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2017_B.txt
#python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2017_C.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2017_D.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2017_E.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2017_F.txt

# 2018
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2018_A.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2018_B.txt
python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2018_C.txt
#python scripts/submit_batch.py -N -1 -n 1 --runtime 6 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/ data/CharmoniumUL_2018_D.txt


