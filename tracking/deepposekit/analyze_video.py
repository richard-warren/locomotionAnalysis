import numpy as np
from deepposekit.models import load_model
from deepposekit.io import DataGenerator, VideoReader, VideoWriter
import pandas as pd
import os

'''
todo:
user input model name, skeleton name, predictions file name, run/wisk
'''

# settings
# max_frames = 50000
video_name = r'Y:\obstacleData\sessions\191110_004\run_short.mp4'
model_name = r'D:\github\obstacle_analysis\tracking\deepposekit\models\model_run.h5'
skeleton = r'D:\github\obstacle_analysis\tracking\label\training_sets\skeleton_run.csv'
output_name = 'trackedFeatures_run.csv'
batch_size = 32

# load model and video
model = load_model(model_name)
reader = VideoReader(video_name, batch_size=batch_size, gray=False, pad_imgs=True)

# predict
# max_batches = max_frames//batch_size
predictions = model.predict(reader, verbose=1)
reader.close()

# save data
file_name = os.path.join(os.path.split(video_name)[0], output_name)
features = list(pd.read_csv(skeleton).name)
columns = np.repeat(features, 3)
data = pd.DataFrame(columns=columns, index=np.arange(predictions.shape[0]))
data[:] = np.reshape(predictions, (predictions.shape[0], -1))
data.to_csv(file_name)