import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

models = [r'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5',
          r'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_DeepLabCut.h5']

# initializations
plt.close()
fig, axes = plt.subplots(1, len(models), figsize=(9,4))
paw_pairs = {'lh': [0,8],
             'lf': [1,9],
             'rf': [2,10],
             'rh': [3,11]}

for i, model in enumerate(models):
    log = os.path.splitext(model)[0] + '_log.h5'

    # get true and predicted values
    with h5py.File(log, 'r') as f:
        y_pred = f['logs']['y_pred'][()][1]
        y_error = f['logs']['y_error'][()][1]
        y = y_error + y_pred  # (image num) X (feature) X (xy)
        euclidean = f['logs']['euclidean'][()][1]

        # make sure the extreme values don't screw us up
        y[y < 0] = np.nan
        euclidean[euclidean > 1000] = np.nan

    # print summary stats
    euclidean_error = np.nanmean(euclidean, axis=0)  # average across frames, then features
    axes[i].annotate('euclidean error: %.2f' % np.nanmean(euclidean_error),
                     xy=(.05, .95), xycoords='axes fraction')
    euclidean_error_paws = euclidean_error[sum(paw_pairs.values(), [])]
    axes[i].annotate('euclidean error (paws only): %.2f' % np.nanmean(euclidean_error_paws),
                     xy=(.05, .9), xycoords='axes fraction')

    # scatter x values for paw in top and bottom views
    for paw, pair in paw_pairs.items():
        axes[i].scatter(y[:, pair[1], 0], y[:, pair[0], 0])
    axes[i].set_title(os.path.split(model)[-1][:-3])




# # find indices for paws
# skeleton = r'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv'
# s = pd.read_csv(skeleton)
# features = list(s.name)

##



##
