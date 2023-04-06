import numpy as np
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentParser
#from cmsjson import CMSJson
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument(
   '--what'
)
parser.add_argument(
   '--test', action='store_true'
)
parser.add_argument(
   '--jobtag', default='prova', type=str
)
parser.add_argument(
   '--ntrees', default=700, type=int               # chiara: era 100. Io sempre usato 500; Rob 2000; Per modello finale 700
)

parser.add_argument(
   '--depth', default=13, type=int                  # chiara: era 4; def = 6; Io ho sempre usato 10; Rob/Mauro: 15; Per modello finale 13
)
parser.add_argument(
   '--lrate', default=0.1, type=float     
)
parser.add_argument(
   '--rstate', default=42, type=int
)
parser.add_argument(
   '--gamma', default=0., type=float
)
parser.add_argument(
   '--min_child_weight', default=1.0, type=int
)
parser.add_argument(
   '--subsample', default=1., type=float
)
parser.add_argument(
   '--colsample_bytree', default=1.0, type=float
)

parser.add_argument(
   '--reg_alpha', default=0.0, type=float
)
parser.add_argument(
   '--reg_lambda', default=2.112612055963768, type=float          # chiara: default 1, sempre usato nei miei test; mauro 9.99999999862; rob: 2.112612055963768
)
parser.add_argument(
   '--nthreads', default=8, type=int
)
parser.add_argument(
   '--no_early_stop', action='store_true'
)
parser.add_argument(
   '--config'
)
parser.add_argument(
   '--dataset'
)
parser.add_argument(
   '--selection'
)
parser.add_argument(
   '--as_weight'
)
parser.add_argument(
   '--noweight', action='store_true'
)
#parser.add_argument(
#   '--SW94X', action='store_true'
#)
parser.add_argument(
   '--usenomatch', action='store_true'
)
parser.add_argument(
   '--load_model', action='store_true'
)
parser.add_argument(
   '--notraining', action='store_true'
)

args = parser.parse_args()

import json
#if args.config:
#   #config overrides eveything
#   cfg = json.load(open(args.config))
#   args.reg_alpha = cfg['reg_alpha'] 
#   args.colsample_bytree = cfg['colsample_bytree'] 
#   args.lrate = cfg['learning_rate'] 
#   args.min_child_weight = cfg['min_child_weight'] 
#   args.ntrees = cfg['n_estimators'] 
#   args.subsample = cfg['subsample'] 
#   args.reg_lambda = cfg['reg_lambda'] 
#   args.depth = cfg['max_depth'] 
#   args.gamma = cfg['gamma']

import matplotlib.pyplot as plt
import ROOT
import uproot
import rootpy
import pandas as pd
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
from datasets import tag, pre_process_data, target_dataset, get_models_dir, reduce_mem_usage
import os

dataset = 'test' if args.test else target_dataset
if args.dataset:
   dataset = args.dataset

mods = '%s/bdt_%s' % (get_models_dir(), args.what)
if not os.path.isdir(mods):
   os.makedirs(mods)
   print' + created directory ', mods

plots = 'MVA/plots/%s/' % (tag)
if not os.path.isdir(plots):
   os.makedirs(plots)
   print' + created directory ', plots 

#from features import *
#features, additional = get_features(args.what)
spectators = ['event','M_X3872', 'M_B0', 'M_PiPi']
features = ['pTM_B0', 'LxySignBSz_B0', 'SVprob_B0', 'CosAlpha3DBSz_B0', 'LxySignSV_K0s', 'SVprob_PiPi', 'pT_PiPi', 'pT_Pi1', 'DR_B0Pi1', 'D0_Pi1'] 

fields = spectators + features
print(fields)
#if args.SW94X and 'seeding' in args.what:
#   fields += seed_94X_additional
#else:
#   fields += additional

#if 'gsf_pt' not in fields : fields += ['gsf_pt'] #@@ redundant?

if not dataset.endswith('.hdf'): # if not args.load_model :
   data = pre_process_data(
      dataset, fields) 
   print(type(data))
      #for_seeding=('seeding' in args.what),
      #keep_nonmatch=args.usenomatch
      #)

   orig = data.copy()                     # all electrons
   print "orig.shape",orig.shape

   if args.selection:
      data = data.query(args.selection)
      print "data.shape", data.shape

   if args.as_weight:
      data['weight'] = data[args.as_weight]

   if args.noweight:
      data['weight'] = 1

   # chiara   
   ####reduce_mem_usage(data)
   from sklearn.model_selection import train_test_split
   #train_test, validation = train_test_split(data, 10, 8)
   train_test, validation = train_test_split(data, test_size=0.2, random_state = 4) 
   #train, test = train_test_split(train_test, 10, 6)
   train, test = train_test_split(train_test, test_size = 0.4, random_state = 4)
   print 'train/test-set ',train_test.shape
   print 'train-set ',train.shape
   print 'test-set ',test.shape
   print 'validation-set ',validation.shape
   print data.is_signal.mean()


from sklearn.externals import joblib
import xgboost as xgb

#
# Train BDTs
#

clf = None
if args.notraining :
   print 'No training done, no pre-existing model loaded!'
elif not args.load_model :

   print 'Training'
   print 'Input features:\n',features

   clf = xgb.XGBClassifier(
      # general parameters
      booster='gbtree',                                # chiara: preso dalla versione di Rob, non nel master (ma e' il default)
      silent=False,
      #### nthread=args.nthreads,
      # booster parameters
      n_estimators=250,
      learning_rate=0.075,
      min_child_weight=args.min_child_weight,          #def=1 in xgboost, as here
      max_depth=3,
      gamma=args.gamma,                                #def=0 in xgboost, as here 
      max_delta_step=0,                                #def=0 in xgboost, as here 
      subsample=args.subsample,                        #def=1 in xgboost, as here     ===> tizio dice che e' tipico iniziare con 0.8
      colsample_bytree=args.colsample_bytree,          #def=1 in xgboost, as here     ===> tizio dice che e' tipico iniziare con 0.8    
      colsample_bylevel=1,                             #def=1; use subsample and colsample_bytree instead
      reg_lambda=args.reg_lambda,                      #def=1 in xgboost, as here  
      reg_alpha=args.reg_alpha,                        #def=0 in xgboost, as here
      scale_pos_weight=1,                              #def=1 in xgboost, as here  
      # learning task parameters
      objective='binary:logitraw',       ## chiara ##
      random_state=5
   )
   
   early_stop_kwargs = {
      'eval_set' : [(test[features].as_matrix(), test.is_signal.as_matrix().astype(int))],
      #'sample_weight_eval_set' : [test.weight.as_matrix()], #undefined in this version
      'eval_metric' : 'auc',
      'early_stopping_rounds' : 15
   } if not args.no_early_stop else {}

   to_train = xgb.DMatrix(train[features].values, label = train.is_signal.values.astype(int),feature_names=features)
   clf.fit(
      train[features].as_matrix(), 
      train.is_signal.as_matrix().astype(int), 
      #####xgb_model='/tmp/crovelli/models/checkpoints/2020Jun5__cmssw_mva_id_nnclean2_BDT.pkl',
      **early_stop_kwargs
   )

   full_model = '%s/%s_%s_%s_BDT.pkl' % (mods, dataset, args.jobtag, args.what)
   joblib.dump(clf, full_model, compress=True)

   print 'Training done!'

else :
   
   full_model = '%s/%s_%s_%s_BDT.pkl' % (mods, dataset, args.jobtag, args.what)
   clf = joblib.load(full_model)
   print 'Loaded pre-existing model!'

# make predictions on test data
pred = clf.predict(test[features].as_matrix())
# compute and print accuracy score
from sklearn.metrics import accuracy_score
print('XGBoost model accuracy score: {0:0.4f}'. format(accuracy_score(test.is_signal.as_matrix().astype(int), pred)))

#
# plot performance
#
from sklearn.metrics import roc_curve, roc_auc_score
args_dict = args.__dict__
plt.rcParams['axes.facecolor'] = 'white'

rocs = {}
if not args.notraining :
   for df, name in [
      ##(train, 'train'),
      ##(test, 'test'),
      (validation, 'validation')
      ]:
      training_out = clf.predict_proba(df[features].as_matrix())[:, 1]
      df['training_out'] = training_out      # chiara: preso dalla versione di Rob, non nel master
      rocs[name] = roc_curve(
            df.is_signal.as_matrix().astype(int), 
            training_out)[:2]
      args_dict['%s_AUC' % name] = roc_auc_score(df.is_signal, training_out)


train_data_out = clf.predict_proba(train.loc[train.is_signal == 0, features].values)[:, 1]
train_MC_out = clf.predict_proba(train.loc[train.is_signal == 1, features].values)[:, 1]
val_data_out = clf.predict_proba(validation.loc[validation.is_signal == 0, features].values)[:, 1]
val_MC_out = clf.predict_proba(validation.loc[validation.is_signal == 1, features].values)[:, 1]

   #with open('%s/%s_%s_%s_BDT.json' % (mods, dataset, args.jobtag, args.what), 'w') as info:
   #   json.dump(args_dict, info)

# make plots
print "Making plots ..."
plt.figure(figsize=[8, 8])
ax = plt.subplot(111)  
box = ax.get_position()   
ax.set_position([box.x0, box.y0, box.width, box.height]) 

plt.title('%s training' % args.what.replace("_"," "))
plt.plot(
   np.arange(0,1,0.01),
   np.arange(0,1,0.01),
   'k--')

plt.plot(rocs['validation'][0][:-1], rocs['validation'][1][:-1], 
            linestyle='solid', 
            color='red', 
            label='UL 2017,  AUC: %.3f'  % args_dict['validation_AUC'])


plt.xlabel('Mistag Rate')
plt.ylabel('Efficiency')
plt.legend(loc='best')
plt.xlim(0., 1)

try : plt.savefig('%s/%s_%s_%s_BDT.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.gca().set_xscale('log')
plt.xlim(1e-4, 1)
try : plt.savefig('%s/%s_%s_%s_log_BDT.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_log_BDT.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
plt.clf()


## feature correlation 
import seaborn as sns
sns.axes_style("whitegrid")
plt.figure(figsize = [10,8])
cov = data[features].corr()
print cov
sns.heatmap(cov, annot=True, fmt=".1f", cmap = 'viridis')
try : plt.savefig('%s/%s_%s_%s_CovMatrix.png' % (plots, dataset, args.jobtag, args.what))
except : pass


## feature imporatance
fRanking= pd.Series(clf.feature_importances_, index = features ).sort_values()
#print(fRanking)
plt.figure(figsize = [8,8])
fRanking.plot(kind='barh')
plt.xlabel('score')
plt.ylabel('features')
plt.grid()
try : plt.savefig('%s/%s_%s_%s_features.png' % (plots, dataset, args.jobtag, args.what))
except : pass

plt.figure(figsize = [8,8])
ax = plt.subplot(111)
plt.title('training VS validation BDT-output ')
nbins =24; xlow = -8; xhigh = 4 
bins = np.linspace(xlow, xhigh, nbins)
counts_data, bin_edges_data = np.histogram(val_data_out, bins, density = True)
counts_MC, bin_edges_MC = np.histogram(val_MC_out, bins, density = True)
bin_centers_MC = (bin_edges_MC[:-1]+bin_edges_MC[1:])/2.
bin_centers_data = (bin_edges_data[:-1]+bin_edges_data[1:])/2.
print bin_centers_MC
print bin_centers_data
ax.hist(train_data_out, bins, density = True, color = 'red', alpha = 0.4, label = 'training data 2017')
ax.hist(train_MC_out  , bins, density = True, color = 'blue', alpha = 0.4, label = 'training MC 2017')
ax.errorbar(bin_centers_data, counts_data, yerr = np.sqrt(counts_data/len(val_data_out)), fmt = 'o', color = 'red', label = 'validation data 2017')
ax.errorbar(bin_centers_MC, counts_MC, yerr = np.sqrt(counts_MC/len(val_MC_out)), fmt = 'o', color = 'blue', label = 'validation MC 2017')
ax.set_xlim(xlow,xhigh)
ax.legend(fontsize = 14)
ax.set_xlabel('BDT output')
ax.grid()
try : plt.savefig('%s/%s_%s_%s_outcome.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_outcome.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass


## BDT output vs Mpipi
plt.figure(figsize = [8,8])
ax = plt.subplot(111)
nbins_BDT = 48; nbins_Mpipi = 40
data2D = ax.hist2d(val_data_out, validation.loc[validation.is_signal == 0, 'M_PiPi'].values, bins = [nbins_BDT, nbins_Mpipi], range = [[xlow,xhigh], [0.2, 1.]], cmap=plt.cm.Reds, cmin=1, label = 'BKG-data (validation set)')
MC2D = ax.hist2d(val_MC_out, validation.loc[validation.is_signal == 1, 'M_PiPi'].values, bins = [nbins_BDT, nbins_Mpipi], range = [[xlow,xhigh], [0.2, 1.]], cmap=plt.cm.Blues, cmin=1, label = 'SGN-MC (validation set)')
ax.set_xlabel("BDT output")
ax.set_ylabel("M(#pi^+ #pi^-)")
ax.legend(['BKG-data (validation set)', 'SGN-MC (validation set)'])
plt.tight_layout()
try : plt.savefig('%s/%s_%s_%s_BDTvsMpipi.png' % (plots, dataset, args.jobtag, args.what))
except : pass
try : plt.savefig('%s/%s_%s_%s_BDTvsMpipi.pdf' % (plots, dataset, args.jobtag, args.what))
except : pass
