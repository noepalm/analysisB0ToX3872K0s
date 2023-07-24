from glob import glob
import os
from pdb import set_trace
#A single place where to bookkeep the dataset file locations
#
tag = ''
posix = ''
target_dataset = ''

import socket
path = 'HLTemulation'

# datasets
input_files = {
    # - 2016 preVFP-
   'UL_2016preVFP_data'  :'/eos/user/c/cbasile/B0toX3872K0s/data/CharmoniumUL_2016preVFP_blind.root',
   'UL_2016preVFP_MC_X'  :'../PRELIMINARY/outRoot/RecoDecay_X3872_UL16preVFP.root',
   'UL_2016preVFP_MC_Psi':'../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16preVFP.root',
    # - 2016 postVFP
   'UL_2016postVFP_data'  :'/eos/user/c/cbasile/B0toX3872K0s/data/CharmoniumUL_2016postVFP_blind.root',
   'UL_2016postVFP_MC_X'  :'../PRELIMINARY/outRoot/RecoDecay_X3872_UL16.root',
   'UL_2016postVFP_MC_Psi':'../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16.root',
    # - 2017 -
   'UL_2017_data'  :'/eos/user/c/cbasile/B0toX3872K0s/data/CharmoniumUL_17_HLTemulation_blind.root',
   'UL_2017_MC_X'  :'../PRELIMINARY/outRoot/RecoDecay_X3872_UL17.root',
   'UL_2017_MC_Psi':'../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL17.root',
    # - 2018 -
   'UL_2018_data'  :'/eos/user/c/cbasile/B0toX3872K0s/data/CharmoniumUL_2018_blind.root',
   'UL_2018_MC_X'  :'../PRELIMINARY/outRoot/RecoDecay_X3872_UL18.root',
   'UL_2018_MC_Psi':'../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL18.root'
}

#import concurrent.futures
import multiprocessing
import uproot
import numpy as np

def get_models_dir():
   
   mods = 'MVA/models/%s/' % (tag)
   if not os.path.isdir(mods):
      os.makedirs(mods)
   return mods

def get_data_sync(dataset, columns, nthreads=2*multiprocessing.cpu_count(), exclude={}, path='HLTemulation'):
   if dataset not in input_files:
      raise ValueError('The dataset %s does not exist, I have %s' % (dataset, ', '.join(input_files.keys())))
   print ' + getting files from "%s": ' % input_files[dataset]
   infiles = uproot.open(input_files[dataset])
   print '  available branches:\n',infiles[path].keys()
   if columns == 'all':
      columns = [i for i in infiles[path].keys() if i not in exclude]
   try:
      ret = infiles[path].arrays(columns)
   except KeyError as ex:
      print 'Exception! ', ex
      set_trace()
      raise RuntimeError("Failed to open %s properly" %(input_files[dataset]))
   return ret

import ROOT
def PDtoROOT_DataFrame(df, columns):
   data = {key: df[key].values for key in columns}
   print(data)
   rdf = ROOT.RDF.FromNumpy(data) # NON funziona grr
   
   return rdf


from sklearn.cluster import KMeans
from sklearn.externals import joblib
import json
from pdb import set_trace

#apply_weight = np.vectorize(lambda x, y: y.get(x), excluded={2})

def kmeans_weighter(features, fname):
   kmeans = joblib.load(fname)
   cluster = kmeans.predict(features)
   str_weights = json.load(open(fname.replace('.pkl', '.json')))   
   weights = {}
   for i in str_weights:
      try:
         weights[int(i)] = str_weights[i]
      except:
         pass
   return apply_weight(cluster, weights)

def training_selection(df,low=0.5,high=15.):
   #'ensures there is a GSF Track and a KTF track within eta/pt boundaries'
   return (df.gsf_mode_pt > low) & (np.abs(df.gsf_mode_eta) < 2.4) & ( (df.gen_dR<=0.03) | (df.gen_dR>=0.1) ) 

import rootpy.plotting as rplt
import root_numpy

class HistWeighter(object):
   def __init__(self, fname):
      values = [[float(i) for i in j.split()] for j in open(fname)]
      vals = np.array(values)
      xs = sorted(list(set(vals[:,0])))
      ys = sorted(list(set(vals[:,1])))
      vals[:,0] += 0.0001
      vals[:,1] += 0.0001
      mask = (vals[:,2] == 0)
      vals[:,2][mask] = 1 #turn zeros into ones
      vals[:,2] = 1/vals[:,2]
      self._hist = rplt.Hist2D(xs, ys)
      root_numpy.fill_hist(self._hist, vals[:,[0,1]], vals[:, 2])

   def _get_weight(self, x, y):
      ix = self._hist.xaxis.FindFixBin(x)
      iy = self._hist.yaxis.FindFixBin(y)
      return self._hist.GetBinContent(ix, iy)

   def get_weight(self, x, y):
      cnt = lambda x, y: self._get_weight(x, y)
      cnt = np.vectorize(cnt)
      return cnt(x, y)

import pandas as pd
import numpy as np
def pre_process_data(dataset, features, training = True, channel = 'X3872', keep_nonmatch=False):  
   
   selection_JpsiPiPi = '(M_X3872 < 3.65831) or (M_X3872 > 3.71493 and M_X3872 < 3.82562) or (M_X3872 > 3.91839)' 
   selection_B0 = '(M_B0 > 5.09783 and M_B0 < 5.189) or (M_B0 > 5.371 and M_B0 < 5.463)'
   mods = get_models_dir()
   
   ## get the BACKGROUND data 
   bkg_dataset = dataset + '_data'
   data_dict = get_data_sync(bkg_dataset, features) 
   bkg_data = pd.DataFrame(data_dict)
   bkg_data['is_signal'] = np.zeros(bkg_data.event.shape).astype(int)
   # ... apply blind region cuts
   if training :
      bkg_data = bkg_data.query(selection_B0)
      bkg_data = bkg_data.query(selection_JpsiPiPi)

   ## get the SIGNAL data 
   if (channel == "X3872"): sgn_dataset = dataset + '_MC_X'
   if (channel == "Psi2S"): sgn_dataset = dataset + '_MC_Psi'
   data_dict = get_data_sync(sgn_dataset, features)
   sgn_data = pd.DataFrame(data_dict)
   sgn_data['is_signal'] = np.ones(sgn_data.event.shape).astype(int)

   ## concatenate and shuffle the two sets
   data = pd.concat([bkg_data, sgn_data])
   # shuffle the DataFrame rows
   data = data.sample(frac = 1, random_state = 33).reset_index(drop=True)
   # print(data.head())
   # print(data.tail())

   return data

def train_analysis_split(data, Ncv):
   mask = data.index % Ncv
   CVsplit = []
   [CVsplit.append(data[mask == i]) for i in range(Ncv)]
   return CVsplit

def reduce_mem_usage(df):
    """ 
    iterate through all the columns of a dataframe and 
    modify the data type to reduce memory usage.        
    """
    start_mem = df.memory_usage().sum() / 1024**2
    print(('Memory usage of dataframe is {:.2f}' 
                     'MB').format(start_mem))

    print 'before'
    print df
    
    for col in df.columns:
        col_type = df[col].dtype
        
        print col

        if col_type != object:
            c_min = df[col].min()
            c_max = df[col].max()
            if str(col_type)[:3] == 'int':
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)  
            else:
                if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                    df[col] = df[col].astype(np.float16)
                elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
        else:
            df[col] = df[col].astype('category')    

    end_mem = df.memory_usage().sum() / 1024**2
    print(('Memory usage after optimization is: {:.2f}' 
                              'MB').format(end_mem))
    print('Decreased by {:.1f}%'.format(100 * (start_mem - end_mem) 
                                             / start_mem))
    
    print 'after'
    print df

    return df
