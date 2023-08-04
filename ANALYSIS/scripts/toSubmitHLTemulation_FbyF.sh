##
#  script to submit jobs to perform the HLT-emulation on data file by file
##

# 2016 PRE VFP
#python scripts/submit_batch.py -N -1 -n 1 -F 0 --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_Bv1.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 2691 --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_Bv2.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 888  --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_C.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 1470 --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_D.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 1246 --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_E.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 785  --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016preVFP_F.txt

# 2016 POST VFP
#python scripts/submit_batch.py -N -1 -n 1 -F 126  --runtime 2 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_F.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 2155 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2016_G.txt

# 2017
python scripts/submit_batch_FilebyFile.py -N -1 -n 1 -F 1255 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/2017B/  -i data/CharmoniumUL_2017_B.txt
python scripts/submit_batch_FilebyFile.py -N -1 -n 1 -F 2634 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/2017C/  -i data/CharmoniumUL_2017_C.txt
python scripts/submit_batch_FilebyFile.py -N -1 -n 1 -F 1369 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/2017D/  -i data/CharmoniumUL_2017_D.txt
python scripts/submit_batch_FilebyFile.py -N -1 -n 1 -F 2176 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/2017E/  -i data/CharmoniumUL_2017_E.txt
python scripts/submit_batch_FilebyFile.py -N -1 -n 1 -F 2945 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/2017F/  -i data/CharmoniumUL_2017_F.txt

# 2018
#python scripts/submit_batch.py -N -1 -n 1 -F 2821 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2018_A.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 1116 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2018_B.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 1314 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2018_C.txt
#python scripts/submit_batch.py -N -1 -n 1 -F 6425 --runtime 12 --eos /eos/user/c/cbasile/B0toX3872K0s/data/jobs/  data/CharmoniumUL_2018_D.txt


