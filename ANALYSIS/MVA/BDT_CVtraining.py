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
   '--jobtag', default='CVtraining', type=str
)
parser.add_argument(
   '--CV', default=3, type=int
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
   '--reg_lambda', default=1., type=float          # chiara: default 1, sempre usato nei miei test; mauro 9.99999999862; rob: 2.112612055963768
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
   '--channel', default = 'X3872', type = str
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
import matplotlib.pyplot as plt
import ROOT
import uproot
import rootpy
import pandas as pd
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
from datasets import tag, pre_process_data, target_dataset, get_models_dir, train_analysis_split, PDtoROOT_DataFrame
import os

if args.dataset:
   dataset = args.dataset

mods = '%s/BDTxgb_%s' % (get_models_dir(), dataset)
if not os.path.isdir(mods):
   os.makedirs(mods)
   print' + created directory ', mods

plots = '/eos/user/c/cbasile/www/B0toX3872K0s/MVA/CVtraining_%s/%s/' % (dataset, tag)
if not os.path.isdir(plots):
   os.makedirs(plots)
   print' + created directory ', plots 

results = 'MVA/results'
if not os.path.isdir(results):
   os.makedirs(results)
   print' + created directory ', results 

spectators = ['event','M_X3872', 'M_B0', 'M_PiPi', 'M_K0s']
features = ['pTM_B0', 'LxySignBSz_B0', 'SVprob_B0', 'CosAlpha3DBSz_B0', 'LxySignSV_K0s', 'SVprob_PiPi', 'pT_PiPi', 'DR_B0Pi1', 'D0_Pi1'] 

fields = spectators + features
print(fields)


# load data from TTree
training_cuts = not args.load_model
if training_cuts: print(" [TRAINING CUTS on data]")
else: print(" [NO CUTS on data]")
data = pre_process_data(dataset, fields, training_cuts, args.channel) 
print(type(data))
   
orig = data.copy()                    
print "orig.shape",orig.shape

if args.selection:
   data = data.query(args.selection)
   print "data.shape", data.shape

if args.as_weight:
   data['weight'] = data[args.as_weight]

if args.noweight:
   data['weight'] = 1

   
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from sklearn.metrics import accuracy_score, roc_curve, roc_auc_score
import xgboost as xgb


#####    BDT structure   #####
clf = None
clf = xgb.XGBClassifier(
   # general parameters
   booster='gbtree',                                # chiara: preso dalla versione di Rob, non nel master (ma e' il default)
   silent=False,
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


#####    CV split  #####
Ncv = args.CV
dataCV = train_analysis_split(data, Ncv)

# -> results
fRanking = []
accuracyTest = []
aucTest = []
rocTest = []
aucScore = []
BDTout = []
rocAnalysis = []
outAnaliysis_df = []
BDToutTrain = []
BDToutTest = []
outTrain_df = []

for i in range(Ncv):
   # Train BDTs on Ncv-1 batch 
   tmp = list(dataCV)
   print(" - element in data %d"%(len(tmp)))
   to_apply = tmp.pop(i)
   print(" - element in data %d"%(len(tmp)))
   train_test = pd.concat(tmp)
   train, test = train_test_split(train_test, test_size = 0.25, random_state = 4)
   print 'train/test-set ',train_test.shape
   print 'train-set ',train.shape
   print 'test-set ',test.shape
   print 'to_apply-set ',to_apply.shape

   if not args.load_model :

      #
      # Train BDTs
      #

      print 'Training set number %d'%(i)
      print 'Input features:\n',features

      early_stop_kwargs = {
         'eval_set' : [(test[features].values, test.is_signal.values)],
         'eval_metric' : 'auc',
         'early_stopping_rounds' : 15
      } if not args.no_early_stop else {}

      clf.fit(
         train[features].values,
         train.is_signal.values, 
         **early_stop_kwargs
      )

      ## save model
      full_model = '%s/%s_BDTtraining_CV%d-%d.pkl' % (mods, dataset, Ncv, i)
      joblib.dump(clf, full_model, compress=True)

      print 'Training done!'


      ## predictions on TRAINING_SET
      train['BDTout'] = clf.predict_proba(train[features].values)[:,1]
      outTrain_df.append(train.copy())
      ## predictions on TEST-SET
      BDToutTest.append(clf.predict_proba(test[features].values)[:,1])
      accuracyTest.append(accuracy_score(test.is_signal.values, 
                                       clf.predict(test[features].values)))
      aucTest.append(roc_auc_score(test.is_signal, BDToutTest[i]))
      rocTest.append(
         roc_curve(
               test.is_signal.values, 
               BDToutTest[i])[:2]
      )
      fRanking.append(pd.Series(clf.feature_importances_, index = features ).sort_values())
      print('XGBoost model accuracy score (test-set): {0:0.4f}'. format(accuracyTest[i]))
      print('XGBoost model    AUC   score (test-set): {0:0.4f}'. format(aucTest[i]))

   else:

      full_model = '%s/%s_BDTtraining_CV%d-%d.pkl' % (mods, dataset, Ncv, i)
      clf = joblib.load(full_model)
      print 'Loaded pre-existing model %s' %(full_model)

   ## predictions on ANALYSIS-SET
   BDTout.append(clf.predict_proba(to_apply[features].values)[:,1])
   data.loc[to_apply.index.values,'BDTout'] = clf.predict_proba(to_apply[features].values)[:,1]
   to_apply['BDTout'] = clf.predict_proba(to_apply[features].values)[:,1]
   outAnaliysis_df.append(to_apply.copy())
   rocAnalysis.append(
      roc_curve(
            to_apply.is_signal.values, 
            BDTout[i])[:2]
   )
   aucScore.append(roc_auc_score(to_apply.is_signal, BDTout[i]))
   print('XGBoost model    AUC   score (analysis-set): {0:0.4f}'. format(aucScore[i]))



#print(data)
import ROOT
## --> SAVE DATA FRAME IN A ROOT TREE
outPath = '%s/%s_BDTtraining_CV%d_%s.csv' %(results, dataset, Ncv, args.channel)
if args.load_model : outPath = '%s/%s_BDTapplication_CV%d_%s.csv' %(results, dataset, Ncv, args.channel) 
data.to_csv(outPath)
print(" [OUT] save pd.DataFrame in %s" %(outPath))
rdf = ROOT.RDF.MakeCsvDataFrame(outPath)
outPath = outPath.replace("csv", "root")
rdf.Snapshot("CVtrainig_%s"%(dataset), outPath)
print(" [OUT] save ROOT-DataFrame in %s" %(outPath))


# # #     # # #
#    PLOTS    #
# # #     # # #

if not args.load_model:


   Nside_x = int(np.sqrt(Ncv))+1
   Nside_y = 2

   ## FEATURES IMPORTANCE ##
   fig, axs = plt.subplots(Nside_x, Nside_y, figsize = [8*Nside_x,6*Nside_y])
   k=0
   for i,j in np.ndindex((Nside_x,Nside_y)):
      if k == Ncv:break
      fRanking[k].plot(kind='barh', ax = axs[i,j])
      k+=1
      #axs[i,j].plot(fRanking[i])
      axs[i,j].set_xlabel('score')
      axs[i,j].set_ylabel('features')
      axs[i,j].grid()
   try : plt.savefig('%s/%s_%s_%s_features.png' % (plots, dataset, args.jobtag, args.what))
   except : pass
   try : plt.savefig('%s/%s_%s_%s_features.pdf' % (plots, dataset, args.jobtag, args.what))
   except : pass

   ## ROC curve for TEST-set ##
   fig, ax = plt.subplots(figsize = [6,6])
   [ax.plot(rocTest[i][0][:-1], rocTest[i][1][:-1], 
               linestyle='solid', 
               #color='red', 
               label='UL-2017 CV %d  AUC: %.3f'  %(i, aucTest[i]))

   for i in range(Ncv)]
   plt.legend()
   plt.grid()
   try : plt.savefig('%s/%s_%s_%s_AUCtest.png' % (plots, dataset, args.jobtag, args.what))
   except : pass
   try : plt.savefig('%s/%s_%s_%s_AUCtest.pdf' % (plots, dataset, args.jobtag, args.what))
   except : pass

   ## ROC curve for ANALYSIS-set ##
   fig, ax = plt.subplots(figsize = [6,6])
   [ax.plot(rocAnalysis[i][0][:-1], rocAnalysis[i][1][:-1], 
               linestyle='solid', 
               #color='red', 
               label='UL-2017 CV %d  AUC: %.3f'  %(i, aucScore[i]))

   for i in range(Ncv)]
   plt.legend()
   plt.grid()
   try : plt.savefig('%s/%s_%s_%s_AUCanalysis.png' % (plots, dataset, args.jobtag, args.what))
   except : pass
   try : plt.savefig('%s/%s_%s_%s_AUCanalysis.pdf' % (plots, dataset, args.jobtag, args.what))
   except : pass

   ## BDTout TRAINING vs ANALYSIS SET
   nbins =24; xlow = -8; xhigh = 4 
   bins = np.linspace(xlow, xhigh, nbins)

   fig, ax = plt.subplots(Nside_x, Nside_y, figsize = [6*Nside_x,6*Nside_y])
   k=0
   for i,j in np.ndindex((Nside_x,Nside_y)):
   # predictions
      if k == Ncv: break
      train_data_out = outTrain_df[k].loc[outTrain_df[k].is_signal == 0, 'BDTout'].values
      train_MC_out   = outTrain_df[k].loc[outTrain_df[k].is_signal == 1, 'BDTout'].values
      val_data_out   = outAnaliysis_df[k].loc[outAnaliysis_df[k].is_signal == 0, 'BDTout'].values
      val_MC_out     = outAnaliysis_df[k].loc[outAnaliysis_df[k].is_signal == 1, 'BDTout'].values
      k+=1

      counts_data, bin_edges_data = np.histogram(val_data_out, bins, density = True)
      counts_MC, bin_edges_MC = np.histogram(val_MC_out, bins, density = True)
      bin_centers_MC = (bin_edges_MC[:-1]+bin_edges_MC[1:])/2.
      bin_centers_data = (bin_edges_data[:-1]+bin_edges_data[1:])/2.
      ax[i,j].hist(train_data_out, 
                  bins, density = True, color = 'red', alpha = 0.4, label = 'training data 2017')
      ax[i,j].hist(train_MC_out, 
                  bins, density = True, color = 'blue', alpha = 0.4, label = 'training MC 2017')
      ax[i,j].errorbar(bin_centers_data, counts_data, 
                     yerr = np.sqrt(counts_data/len(val_data_out)), 
                     fmt = 'o', color = 'red', label = 'validation data 2017')
      ax[i,j].errorbar(bin_centers_MC, counts_MC, 
                     yerr = np.sqrt(counts_MC/len(val_MC_out)), 
                     fmt = 'o', color = 'blue', label = 'validation MC 2017')
      ax[i,j].set_xlim(xlow,xhigh)
      ax[0,0].legend(fontsize = 14)
      ax[i,j].set_xlabel('BDT output')
      ax[i,j].grid()
   try : plt.savefig('%s/%s_%s_%s_BDTout_trainVSanalysis.png' % (plots, dataset, args.jobtag, args.what))
   except : pass
   try : plt.savefig('%s/%s_%s_%s_BDTout_trainVSanalysis.pdf' % (plots, dataset, args.jobtag, args.what))
   except : pass

   ## BDT output vs Mpipi
   plt.figure(figsize = [8,8])
   ax = plt.subplot(111)
   nbins_BDT = 48; nbins_Mpipi = 40
   data2D = ax.hist2d(data.loc[data.is_signal == 0, 'BDTout'].values, 
                        data.loc[data.is_signal == 0, 'M_PiPi'].values, 
                        bins = [nbins_BDT, nbins_Mpipi], 
                        range = [[xlow,xhigh], [0.2, 1.]], 
                        cmap=plt.cm.Reds, 
                        cmin=1, 
                        label = 'BKG-data (validation set)')
   MC2D   = ax.hist2d(data.loc[data.is_signal == 1, 'BDTout'].values, 
                     data.loc[data.is_signal == 1, 'M_PiPi'].values, 
                     bins = [nbins_BDT, nbins_Mpipi], 
                     range = [[xlow,xhigh], [0.2, 1.]], 
                     cmap=plt.cm.Blues, alpha = 0.8,
                     cmin=1, 
                     label = 'SGN-MC (validation set)')
   ax.set_xlabel("BDT output")
   ax.set_ylabel("M(#pi^+ #pi^-)")
   ax.legend(['BKG-data (validation set)', 'SGN-MC (validation set)'])
   plt.tight_layout()
   try : plt.savefig('%s/%s_%s_%s_BDTvsMpipi.png' % (plots, dataset, args.jobtag, args.what))
   except : pass
   try : plt.savefig('%s/%s_%s_%s_BDTvsMpipi.pdf' % (plots, dataset, args.jobtag, args.what))
   except : pass
