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
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from datasets import tag, pre_process_data, target_dataset, get_models_dir, train_test_split, reduce_mem_usage
import os
#from bayes_opt import BayesianOptimization



parser = ArgumentParser()
parser.add_argument(
    '--what'
)

parser.add_argument(
   '--jobtag', default='prova', type=str
)

parser.add_argument(
   '--dataset' 
)

parser.add_argument(
   '--selection' 
)

parser.add_argument(
   '--test', action='store_true'
)

parser.add_argument(
   '--noweight', action='store_true'
)

parser.add_argument(
   '--recover', action='store_true',      # chiara
   help='recover lost best iteration due to bug'
)


args = parser.parse_args()

# + Load dataset
if args.dataset:
   dataset = args.dataset

mods = '%s/bdt_%s' % (get_models_dir(), args.what)
if not os.path.isdir(mods):
   os.makedirs(mods)
   print( " + created directory " + mods)

plots = 'MVA/plots/%s/' % (tag)
if not os.path.isdir(plots):
   os.makedirs(plots)
   print(' + created directory ', plots) 

additional = ['event','M_X3872', 'M_B0']
features   = ['pTM_B0', 'LxySignBSz_B0', 'SVprob_B0', 'CosAlpha3DBSz_B0', 'LxySignSV_K0s', 'SVprob_PiPi', 'pT_PiPi', 'pT_Pi1', 'DR_B0Pi1', 'D0_Pi1'] 

fields = additional + features

data = pre_process_data(dataset, fields) 
orig = data.copy() 
print(" orig.shape",orig.shape)

## CONTROL PLOTS ##
signal     = data.query("is_signal == 1")
background = data.query("is_signal == 0")

fig, ax= plt.subplots(figsize=(6, 6))
h = ax.hist2d(background.M_B0.values, background.M_X3872.values, bins = (70,50), cmap=plt.cm.jet, cmin=1)
fig.colorbar(h[3], ax=ax)
ax.set_xlabel("M(B0)(GeV)")
ax.set_ylabel("M(JpsiPiPi) (GeV)")
plt.savefig("/eos/user/c/cbasile/www/B0toX3872K0s/MVA/prove/BkgData.png")
plt.savefig("/eos/user/c/cbasile/www/B0toX3872K0s/MVA/prove/BkgData.pdf")
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
hyperparameter_space = {
    'max_depth': range(3, 15, 2),
    'learning_rate': np.arange(0.05, 0.3, 0.025),
    'n_estimators' : range(50,300,50),
}
print(" ### OPTIMIZATION ### ")
print hyperparameter_space
seed = 43
Niter = 50

classifier = xgb.XGBClassifier(objective='binary:logitraw');
optimizer = RandomizedSearchCV(classifier, 
                        param_distributions= hyperparameter_space, 
                        scoring='roc_auc',
                        n_iter = Niter, 
                        cv = 3, 
                        verbose = 3)
#optimizer = GridSearchCV(classifier,
#                param_grid = hyperparameter_space,
#                scoring = 'roc_auc',
#                cv = 3,
#                verbose = 3)
print(" ... starting random search ")
optimizer.fit(train_test[features].values, train_test.is_signal.values.astype(int)) 

print("... random search is over")
print(" = best model is ", optimizer.best_params_)
print(" = AUC score of the best model is %.3f"%optimizer.score(validation[features].values, validation.is_signal.values.astype(int)))

opt_results = pd.DataFrame(optimizer.cv_results_)
fig = plt.figure()
#plt.errorbar(range(0,Niter),opt_results['mean_test_score'], yerr =opt_results['std_test_score'], marker = 'o', mfc = 'red', label='Random search')
plt.plot(opt_results['mean_test_score'], color = 'red', label='Random search')
#plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel('AUC')
plt.grid(True)
plt.legend()
plt.savefig("/eos/user/c/cbasile/www/B0toX3872K0s/MVA/prove/RandomSearch.png")
plt.savefig("/eos/user/c/cbasile/www/B0toX3872K0s/MVA/prove/RandomSearch.pdf")

