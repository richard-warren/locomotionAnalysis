import sys
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import glob
from deepposekit.io import TrainingGenerator, DataGenerator
from deepposekit.augment import FlipAxis
import imgaug.augmenters as iaa
import imgaug as ia
from deepposekit.models import StackedDenseNet, DeepLabCut, StackedHourglass, LEAP
from deepposekit.models import load_model
from tensorflow.keras.callbacks import ReduceLROnPlateau, EarlyStopping
from deepposekit.callbacks import Logger, ModelCheckpoint
import h5py
import time
from os.path import expanduser
import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'

# pass 'run' or 'wisk' as the only command line argument
view = 'run' if len(sys.argv)==1 else sys.argv[1]
print('training %s model' % view)

## make data generator

training_set = r'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_%s.h5' % view
data_generator = DataGenerator(training_set, zeros_to_nan=True)

print('data_generator of length %i and size (%i, %i) created' %
      (len(data_generator), data_generator[0][0].shape[1], data_generator[0][0].shape[2]))

# ## visualize frame
#
# image, keypoints = data_generator[0]
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

## set up data augmentation

augmenter = []
sometimes = []

sometimes.append(iaa.Affine(scale={"x": (0.98, 1.02), "y": (0.98, 1.02)},
                            translate_percent={'x': (-0.1, 0.1), 'y': (-0.1, 0.1)},
                            shear=(-8, 8),
                            order=ia.ALL,
                            cval=ia.ALL,
                            mode=ia.ALL)
                 )
sometimes.append(iaa.Affine(scale=(0.6, 1.4),
                            mode=ia.ALL,
                            order=ia.ALL,
                            cval=ia.ALL)
                 )
augmenter.append(iaa.Sometimes(0.75, sometimes))
augmenter.append(iaa.Affine(rotate=(-5, 5),
                            mode=ia.ALL,
                            order=ia.ALL,
                            cval=ia.ALL)
                 )
augmenter = iaa.Sequential(augmenter)


# ## check out augmentations
# plt.close('all')
# image, keypoints = data_generator[0]
# image, keypoints = augmenter(images=image, keypoints=keypoints)
# plt.figure(figsize=(5,5))
# image = image[0] if image.shape[-1] is 3 else image[0, ..., 0]
# cmap = None if image.shape[-1] is 3 else 'gray'
# plt.imshow(image, cmap=cmap, interpolation='none')
# for idx, jdx in enumerate(data_generator.graph):
#     if jdx > -1:
#         plt.plot(
#             [keypoints[0, idx, 0
#              ], keypoints[0, jdx, 1]],
#             'r-'
#         )
# plt.scatter(keypoints[0, :, 0], keypoints[0, :, 1], c=np.arange(data_generator.keypoints_shape[0]), s=50, cmap=plt.cm.hsv, zorder=3)
# plt.xlim(0, image.shape[1])
# plt.ylim(image.shape[0], 0)
# plt.show()

## make training generator


train_generator = TrainingGenerator(generator=data_generator,
                                    downsample_factor=3,
                                    augmenter=augmenter,
                                    sigma=10,
                                    validation_split=0.2,
                                    use_graph=True,
                                    random_seed=1,
                                    graph_scale=1)
train_generator.get_config()


# ## check training generator output
#
# n_keypoints = data_generator.keypoints_shape[0]
# batch = train_generator(batch_size=1, validation=False)[0]
# inputs = batch[0]
# outputs = batch[1]
#
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
# ax1.set_title('image')
# ax1.imshow(inputs[0,...,0], cmap='gray', vmin=0, vmax=255)
#
# ax2.set_title('posture graph')
# ax2.imshow(outputs[0,...,n_keypoints:-1].max(-1))
#
# ax3.set_title('keypoints confidence')
# ax3.imshow(outputs[0,...,:n_keypoints].max(-1))
#
# ax4.set_title('posture graph and keypoints confidence')
# ax4.imshow(outputs[0,...,-1], vmin=0)
# plt.show()
#
# train_generator.on_epoch_end()

## define model

# model = StackedDenseNet(train_generator, n_stacks=2, growth_rate=32, pretrained=False)
model = DeepLabCut(train_generator, backbone="resnet50")
model.get_config()


# ## test prediction speed
#
# batch_size = 32
#
# data_size = (10000,) + data_generator.image_shape
# x = np.random.randint(0, 255, data_size, dtype="uint8")
# y = model.predict(x[:100], batch_size=batch_size)  # make sure the model is in GPU memory
# t0 = time.time()
# y = model.predict(x, batch_size=batch_size, verbose=1)
# t1 = time.time()
# print('frames per second: %.1f' % (x.shape[0] / (t1 - t0)))

##  set up training

model_folder = os.path.join('tracking', 'deepposekit', 'models')
model_name = 'model_%s_%s' % (view, model.get_config()['name'])

logger = Logger(validation_batch_size=8,
    filepath = os.path.join(model_folder, model_name+'_log.h5')
)

reduce_lr = ReduceLROnPlateau(monitor="val_loss", factor=0.2, verbose=1, patience=20)

model_checkpoint = ModelCheckpoint(
    os.path.join(model_folder, model_name+'.h5'),
    monitor="val_loss",
    # monitor="loss" # use if validation_split=0
    verbose=1,
    save_best_only=True,
)

early_stop = EarlyStopping(
    monitor="val_loss",
    # monitor="loss" # use if validation_split=0
    min_delta=0.001,
    patience=100,
    verbose=1
)

callbacks = [logger, reduce_lr, model_checkpoint, early_stop]
# callbacks = [reduce_lr, model_checkpoint, early_stop]  # !!! does not work with logger for some reason...

## train model

history = model.fit(
    batch_size=8,
    validation_batch_size=8,
    callbacks=callbacks,
    epochs=250,
    n_workers=8,
    steps_per_epoch=None,
)

# plot training history
plt.figure()
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.legend(['Train', 'Test'], loc='upper left')
plt.show()