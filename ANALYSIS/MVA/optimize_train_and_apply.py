import numpy as np
import os
import json
from argparse import ArgumentParser
#from joblib import Parallel, delayed

parser = ArgumentParser()
parser.add_argument(
    '--no_training', action='store_true', default = False 
)

args = parser.parse_args()

# TODO: can be easily parallelized (each era iteration in the for loop is independent)

# def optimize_train_and_apply(era):
    #all code below
#Parallel(n_jobs=4)(delayed(optimize_train_and_apply)(era) for era in ['D', 'E', 'F', 'G'])

for era in ['D', 'E', 'F', 'G']:

    print("####### 2022, ERA %s #######" % era)

    if not args.no_training:
        ##### HYPERPARAMETER OPTIMIZATION #####
        print("### PERFORMING HYPERPARAMETER OPTIMIZATION ON ERA %s ###" % era)

        # Check if hyperparameter optimization has been already performed
        if os.path.exists('MVA/models/BDT_optimisation_random_2022%s/BestBDT_random_2022%s_3.json' % (era, era)):
            print("   RandomSearch for hyperparameter optimization already performed. Moving to GridSearch...")
        else:
            # Run random search 3 times to find grid search range
            for i in range(3):
                print("   RANDOM SEARCH %d" % (i+1))
                os.system("python BDT_HyperParams_optimiser.py --dataset 2022 --era %s" % era)
                os.system("mv MVA/models/BDT_optimisation_random_2022%s/BestBDT_random_2022%s.json MVA/models/BDT_optimisation_random_2022%s/BestBDT_random_2022%s_%d.json" % (era, era, era, era, i+1))

        # Check if grid search has been already performed

        if os.path.exists('MVA/models/BDT_optimisation_grid_2022%s/BestBDT_grid_2022%s.json' % (era, era)):
            print("   GridSearch already performed. Moving to training...")
        else:
            # Retrieve json with best parameter and find range
            learning_rate = []
            n_estimators = []
            max_depth = []

            for i in range(3):
                f = open('MVA/models/BDT_optimisation_random_2022%s/BestBDT_random_2022%s_%d.json' % (era, era, i+1))
                params = json.load(f)
                learning_rate.append(params['learning_rate'])
                n_estimators.append(params['n_estimators'])
                max_depth.append(params['max_depth'])

            print("   Best hyperparameters found per random search:")
            print("     " + json.dumps({'learning_rate': learning_rate, 'n_estimators': n_estimators, 'max_depth': max_depth}))

            n_estimators_gs  = np.arange(min(n_estimators), max(n_estimators)+1, 50).tolist()      #tolist to make it json serializable
            max_depth_gs     = np.arange(min(max_depth), max(max_depth)+1, 1).tolist()
            learning_rate_gs = np.linspace(min(learning_rate), max(learning_rate), 3).tolist()

            # save grid to file to be retrived by hyperparameter optimiser
            grid_search_params = {'learning_rate': learning_rate_gs, 'n_estimators': n_estimators_gs, 'max_depth': max_depth_gs}
            
            ## check if output folder already exists, else create
            grid_outfolder = 'MVA/models/BDT_optimisation_grid_2022%s' % era
            if not os.path.isdir(grid_outfolder):
                os.makedirs(grid_outfolder)
                print("   created directory %s" % grid_outfolder)

            ## dump grid to file
            json.dump(grid_search_params, open('MVA/models/BDT_optimisation_grid_2022%s/BDT_grid_to_search_2022%s.json' % (era, era), 'w'))
            print("   Grid search range:")            
            print("     " + json.dumps(grid_search_params))

            # Perform grid search
            print("   GRID SEARCH")
            os.system("python BDT_HyperParams_optimiser.py --dataset 2022 --era %s --grid_search --from_file" % era)

    ##### TRAINING #####
    
    # Retrieve best parameters
    f = open('MVA/models/BDT_optimisation_grid_2022%s/BestBDT_grid_2022%s.json' % (era, era))
    params = json.load(f)
    print("   [!] Best hyperparameters: %s" % json.dumps(params))
    
    print("### TRAINING AND APPLYING BDT ON ERA %s ###" % era)

    # Train BDT
    # check if model was already trained
    if os.path.exists('/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022%s_BDTtraining_CV3_X3872.root' % era):
        print("   BDT already trained. Moving to application...")
    else:
        print("   TRAINING BDT")
        os.system("python BDT_CVtraining.py --dataset 2022 --era %s --what %s --CV 3 --lrate %f --ntrees %d --depth %d" % (era, era, params['learning_rate'], params['n_estimators'], params['max_depth']))

    # Apply BDT
    ## X3872
    # check if model was already applied
    if os.path.exists('/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022%s_BDTapplication_CV3_X3872.root' % era):
        print("   BDT already applied on X3872. Moving to Psi2S...")
    else:        
        print("   APPLYING BDT ON X3872")
        os.system("python BDT_CVtraining.py --dataset 2022 --era %s --what %s --CV 3 --channel X3872 --load_model" % (era, era))

    ## Psi2S
    # check if model was already applied
    if os.path.exists('/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022%s_BDTapplication_CV3_Psi2S.root' % era):
        print("   BDT already applied on Psi2S. Moving to next era...")
    else:
        print("   APPLYING BDT ON Psi2S")
        os.system("python BDT_CVtraining.py --dataset 2022 --era %s --what %s --CV 3 --channel Psi2S --load_model" % (era, era))

