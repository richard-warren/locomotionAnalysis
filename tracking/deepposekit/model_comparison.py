import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
from sklearn.linear_model import LinearRegression


## COMPARE PERFORMANCE ON VALIDATION SETS

# settings
models = [r'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5',
          r'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_DeepLabCut.h5']
paw_pairs = {'lh': [0,8],
             'lf': [1,9],
             'rf': [2,10],
             'rh': [3,11]}

plt.close()
fig, axes = plt.subplots(1, len(models), figsize=(9,4))

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
        axes[i].scatter(y[:, pair[1], 0], y[:, pair[0], 0], s=8, alpha=.5)
    axes[i].set_title(os.path.split(model)[-1][:-3])

plt.savefig(os.path.join('tracking', 'deepposekit', 'figures', 'errors'))

## COMPARE PERFORMANCE ON FULL VIDEOS

# settings
paw_names = ['paw1LH', 'paw2LF', 'paw3RF', 'paw4RH']  # names of fields in tracking spreadsheets
sessions = ['200113_000', '200116_000', '200117_000', '200114_000', '200131_000', '200202_000', '191221_000']
files = ['trackedFeatures_run.csv', 'trackedFeatures_runDLC.csv']
scatters = 1000
thresh = 25


# make df with one row per tracked feature, for all features in all frames in all sessions
# todo: this is super memory inefficient // perhaps subsample for entire videos... // also, should have different columns for dft sessions, because having to do logical indexing below over this enormous length df takes tons of time
data = pd.DataFrame(columns=['tracking_file', 'session', 'paw', 'x_bot', 'x_top'])

for s in sessions:
    for f in files:
        csv_name = os.path.join(os.environ['OBSDATADIR'], 'sessions', s, f)
        if os.path.exists(csv_name):
            print('%s: loading tracking data...' % s)
            tracking = pd.read_csv(csv_name)

            for paw in paw_names:
                idx = len(data)
                data = data.reindex(list(range(len(data) + len(tracking))))  # extend dataframe

                data.loc[idx:, 'tracking_file'] = f
                data.loc[idx:, 'session'] = s
                data.loc[idx:, 'paw'] = paw
                data.loc[idx:, 'x_bot'] = list(tracking[paw+'_bot'])
                data.loc[idx:, 'x_top'] = list(tracking[paw+'_top'])
        else:
            print('%s: WARNING! %s does not exist for this session...' % (s, f))

print('all done!')

##
plt.close()
fig, axes = plt.subplots(len(files), len(sessions), figsize=(len(sessions)*2.5, len(files)*2.5))
regressor = LinearRegression()
error_rates = np.empty((len(files), len(sessions)))

for s_i, s in enumerate(sessions):
    for f_i, f in enumerate(files):
        inds = (data.session == s) & (data.tracking_file == f)
        for paw in paw_names:
            inds_sub = inds & (data.paw == paw)
            x = np.array(data.loc[inds_sub, 'x_bot']).reshape(-1, 1)
            y = np.array(data.loc[inds_sub, 'x_top'])

            if len(x):
                inds_sub = random.sample(range(len(x)), scatters//len(paw_names))
                axes[f_i, s_i].scatter(x[inds_sub], y[inds_sub], s=.75, alpha=.5)
        axes[f_i, s_i].set_xticks([])
        axes[f_i, s_i].set_yticks([])

        if s_i==0:
            axes[f_i, s_i].set_ylabel(os.path.splitext(f)[0])
        if f_i==0:
            axes[f_i, s_i].set_title(s)

        # find error rate for session
        x = np.array(data.loc[inds, 'x_bot']).reshape(-1, 1)
        y = np.array(data.loc[inds, 'x_top'])
        if len(x):
            regressor.fit(x, y)
            y_hat = regressor.predict(x)
            error_rates[f_i, s_i] = np.mean(np.abs((y - y_hat) > thresh))
            axes[f_i, s_i].annotate('err rate: %.5f' % error_rates[f_i, s_i], xy=(.04, .9), xycoords='axes fraction')

fig.savefig(os.path.join('tracking', 'deepposekit', 'figures', 'x_alignment'))

## print average error rates across models
print('------ERROR RATES------')
for f_i, f in enumerate(files):
    print('%s: %.5f' % (f, np.mean(error_rates[f_i])))

## todo: scatter error rates with two algorithms against eachother (only works if len(files)==2)