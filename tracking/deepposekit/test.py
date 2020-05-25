import numpy as np
import matplotlib.pyplot as plt
from deepposekit.io import TrainingGenerator, DataGenerator

## make data generator

training_set = r'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.h5'
data_generator = DataGenerator(training_set, zeros_to_nan=True)
print('data_generator of length %i created' % len(data_generator))

# find data_generator index where one of the paws is out of frame
for i, (img, k) in enumerate(data_generator):
    if np.any(k[0,:4,0]<0):
        dg_idx = i
        print('data_generator index: %i' % i)
        break


# ## visualize frame
#
# image, keypoints = data_generator[dg_idx]
#
# plt.figure(figsize=(5,5))
# image = image[0] if image.shape[-1] is 3 else image[0, ..., 0]
# cmap = None if image.shape[-1] is 3 else 'gray'
# plt.imshow(image, cmap=cmap, interpolation='none')
# for idx, jdx in enumerate(data_generator.graph):
#     if jdx > -1:
#         plt.plot(
#             [keypoints[0, idx, 0], keypoints[0, jdx, 0]],
#             [keypoints[0, idx, 1], keypoints[0, jdx, 1]],
#             'r-'
#         )
# plt.scatter(keypoints[0, :, 0], keypoints[0, :, 1], c=np.arange(data_generator.keypoints_shape[0]), s=50, cmap=plt.cm.hsv, zorder=3)
# plt.xlim(0, image.shape[1])
# plt.ylim(image.shape[0], 0)
# plt.show()

## make training generator

train_generator = TrainingGenerator(generator = data_generator,
                                    downsample_factor=3,
                                    sigma=10,
                                    validation_split=0,
                                    use_graph=True,
                                    random_seed=1,
                                    graph_scale=1,
                                    shuffle=False)
train_generator.get_config()

## check training generator output
plt.close()
n_keypoints = data_generator.keypoints_shape[0]
batch = train_generator(batch_size=1, validation=False)[dg_idx]
inputs = batch[0]
outputs = batch[1]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
ax1.set_title('image')
ax1.imshow(inputs[0,...,0], cmap='gray', vmin=0, vmax=255)

ax2.set_title('posture graph')
ax2.imshow(outputs[0,...,n_keypoints:-1].max(-1))

ax3.set_title('keypoints confidence')
ax3.imshow(outputs[0,...,:n_keypoints].max(-1))

ax4.set_title('posture graph and keypoints confidence')
ax4.imshow(outputs[0,...,-1], vmin=0)
plt.show()

train_generator.on_epoch_end()

## test data generator cropping

from deepposekit.io import VideoReader
import matplotlib.pyplot as plt
import numpy as np

reader = VideoReader(r'Z:\loco\obstacleData\sessions\200308_000\run_originalDimensions.mp4',
                     batch_size=1, gray=True, frame_size=[600,1000])
plt.imshow(np.squeeze(reader[0]))