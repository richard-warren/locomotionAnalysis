import sys
import numpy as np
from deepposekit.models import load_model
from deepposekit.io import DataGenerator, VideoReader, VideoWriter
import pandas as pd
import os

# parse input arguments
video, model_name, skeleton, output = sys.argv[1:]
# video = r'Z:\loco\obstacleData\sessions\200308_000\run_originalDimensions.mp4'
# model_name = r'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run.h5'
# skeleton = r'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv'
# output = 'trackedFeatures_run.csv'

# settings
batch_size = 32
# max_frames = None  # set to None unless debugging
max_frames = None  # set to None unless debugging

# load model and video
model = load_model(model_name)
reader = VideoReader(video, batch_size=batch_size, gray=True, pad_imgs=True)

# predict
max_batches = max_frames//batch_size if max_frames else None
predictions = model.predict(reader, verbose=1, steps=max_batches)
reader.close()

# save data
file_name = os.path.join(os.path.split(video)[0], output)
features = list(pd.read_csv(skeleton).name)
columns = np.repeat(features, 3)
data = pd.DataFrame(columns=columns, index=np.arange(predictions.shape[0]))
data[:] = np.reshape(predictions, (predictions.shape[0], -1))
data.to_csv(file_name)