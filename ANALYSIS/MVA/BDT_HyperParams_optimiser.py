import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
#from cmsjson import CMSJson
import matplotlib.pyplot as plt
import ROOT
import uproot
import rootpy
import pandas as pd
import json 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import tag, pre_process_data, target_dataset, get_models_dir 
import os
import logging
mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

#######
# -> optimisation of the BDT hyperparameters by random/grid search
#    random-search usage : python BDT_HyperParams_optimiser.py --dataset UL_2017 
#    grid-search usage : python BDT_HyperParams_optimiser.py --dataset UL_2017 --grid_search
#
#    NOEMI: python BDT_HyperParams_optimiser.py --dataset 2022 --era D
#    NOEMI: python BDT_HyperParams_optimiser.py --dataset 2022 --era E
#    NOEMI: python BDT_HyperParams_optimiser.py --dataset 2022 --era F
#    NOEMI: python BDT_HyperParams_optimiser.py --dataset 2022 --era G
#######

parser = ArgumentParser()
#parser.add_argument(
#    '--what'
#)
""
parser.add_argument(
   '--jobtag', default='', type=str
)

parser.add_argument(
   '--dataset' 
)

parser.add_argument(
   '--era' 
)

parser.add_argument(
   '--selection' 
)

parser.add_argument(
   '--grid_search', action='store_true', default = False
)

parser.add_argument(
   '--from_file', action='store_true', default = False
)


#parser.add_argument(
#   '--test', action='store_true'
#)

#parser.add_argument(
#   '--noweight', action='store_true'
#)

#parser.add_argument(
#   '--recover', action='store_true',      # chiara
#   help='recover lost best iteration due to bug'
#)


args = parser.parse_args()

# + Load dataset
if args.dataset:
   dataset = args.dataset

# + Load and check era
if args.era:
    # era must be a single, uppercase letter OR "all"
    if (len(args.era) == 1 and args.era.isupper()) or args.era == "all":
        era = args.era
    else:
        raise ValueError("ERROR: era must be a single, uppercase letter OR 'all'")
else:
    print("ERROR: no era specified")
    exit()

# + Set job tag
if args.grid_search:
    mode = 'grid' 
else :
    mode = 'random'

# mods = '%s/BDT_optimisation_%s_%s%s' % (get_models_dir(), mode, dataset, era)
mods = '%s/BDT_optimisation_%s_%s' % (get_models_dir(), mode, dataset)
if not os.path.isdir(mods):
   os.makedirs(mods)
   print( " + created directory " + mods)

plots = './plots/'
if not os.path.isdir(plots):
   os.makedirs(plots)
   print(' + created plot directory ', plots)

additional = ['event','M_X3872', 'M_B0']
features   = ['pTM_B0', 'LxySignBSz_B0', 'SVprob_B0', 'CosAlpha3DBSz_B0', 'LxySignSV_K0s', 'SVprob_PiPi', 'pT_PiPi', 'pT_Pi1', 'DR_B0Pi1', 'D0_Pi1'] 

fields = additional + features

data = pre_process_data(dataset, fields, era = era, keep_nonmatch = False) # optimize on MC-matched signal only
orig = data.copy() 
print(" orig.shape",orig.shape)

# once data has been processed, change era to "" for output file names
if era == "all": 
      era = ""
print(type(data))

## CONTROL PLOTS ##
signal     = data.query("is_signal == 1")
background = data.query("is_signal == 0")

fig, ax= plt.subplots(figsize=(6, 6))
h = ax.hist2d(background.M_B0.values, background.M_X3872.values, bins = (70,50), cmap=plt.cm.jet, cmin=1)
fig.colorbar(h[3], ax=ax)
ax.set_xlabel("M(B0)(GeV)")
ax.set_ylabel("M(JpsiPiPi) (GeV)")
fig_name = 'MB0vsMJpsiPiPi_%s_%s%s' %(mode, dataset, era)
plt.savefig('%s/%s.png'%(plots, fig_name))
plt.savefig('%s/%s.pdf'%(plots, fig_name))
#exit()

if args.selection:
   data = data.query(args.selection)
   print (" data.shape", data.shape)


from sklearn.model_selection import train_test_split
train_test, validation = train_test_split(data, test_size=0.2)
train, test = train_test_split(train_test, test_size = 0.4)
print (' train/test-set ',train_test.shape)
print (' train-set ',train.shape)
print (' test-set ',test.shape)
print (' validation-set ',validation.shape)
print (data.is_signal.mean())


from sklearn.externals import joblib
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, cross_val_score
from scipy.stats import uniform
import xgboost as xgb

# hyperparameters to be tuned
if not (args.grid_search):
    hyperparameter_space = {
        'max_depth': range(3, 10),
        'learning_rate': np.arange(0.05, 0.3, 0.025),
        'n_estimators' : range(50,500,50),
    }
else :
    if (args.from_file):
        hyperparameter_space = json.load(open('MVA/models/BDT_optimisation_grid_2022%s/BDT_grid_to_search_2022%s.json' % (era, era)))
    else :
        hyperparameter_space = {
            'max_depth': [3,4], 
            'learning_rate': [0.1, 0.125, 0.175], 
            'n_estimators' : range(100,300,50),
        }
    
print(" ### OPTIMIZATION ### ")
print(hyperparameter_space)
seed = 43
Niter = 5 #50

classifier = xgb.XGBClassifier(objective='binary:logitraw');
if not (args.grid_search):
    optimizer = RandomizedSearchCV(classifier, 
                            param_distributions= hyperparameter_space, 
                            scoring='roc_auc',
                            n_iter = Niter, 
                            cv = 3, 
                            verbose = 3,
                            n_jobs = -1) #use all available cores
    print(" ... starting random search ")

if (args.grid_search):
    optimizer = GridSearchCV(classifier,
                    param_grid = hyperparameter_space,
                    scoring = 'roc_auc',
                    cv = 3,
                    verbose = 3,
                    n_jobs = -1) #use all available cores
    print(" ... starting grid search ")

optimizer.fit(train_test[features].values, train_test.is_signal.values.astype(int)) 

print("... hyper parameters search is over")
print(" = best model is ", optimizer.best_params_)
outfile_name = '%s/BestBDT_%s_%s%s.json'%(mods, mode, dataset, era)
with open(outfile_name, 'w') as outfile:
    json.dump(optimizer.best_params_, outfile)
print(" = validation AUC score of the best model is %.3f"%optimizer.best_estimator_.score(validation[features].values, validation.is_signal.values.astype(int)))

opt_results = pd.DataFrame(optimizer.cv_results_)
fig = plt.figure()
#plt.errorbar(range(0,Niter),opt_results['mean_test_score'], yerr =opt_results['std_test_score'], marker = 'o', mfc = 'red', label='Random search')
lab = 'Random search'
if (args.grid_search) : lab = 'Grid search'
plt.plot(opt_results['mean_test_score'], color = 'red', label=lab)
#plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel('AUC')
plt.grid(True)
plt.legend()
fig_name = 'AUCscore_%s_%s%s_%s' %(mode, dataset, era, lab)
plt.savefig('%s/%s.png'%(plots, fig_name))
plt.savefig('%s/%s.pdf'%(plots, fig_name))