from moviepy.editor import VideoFileClip
import keras
import numpy as np
import logging
from skimage.transform import resize
import operator
import keras.backend as K
import scipy.io
import bisect
from keras.models import Model
from keras.layers import Input, MaxPooling2D
import time

import csv

from tqdm import tqdm

from slowai.layers import *
from slowai.metrics import *

import os

import sys

from collections import deque

if len(sys.argv)!=3:
    print("Usage: python cropanalyzevideo.py baseDir session")
    exit()

if os.path.split(sys.argv[0])[0] != '':
    os.chdir(os.path.split(sys.argv[0])[0])

session = sys.argv[2]
baseDir = sys.argv[1]
modelName = "largecrop.20-1.40"
cropModelName = "cropWhiskers"
probDistribution = [0.5, 0.707, 0.867, 0.966, 1.0, 1.0, 0.966, 0.867, 0.707, 0.5]
timesteps = 10
input_size = (280,336)
color_mode = "grayscale"
numFeatures = 5408
bs = 64

print("Loading models")

def meanError(y_true, y_pred):
    groundTruth = K.argmax(y_true, axis=-1)
    preds = K.argmax(y_pred, axis=-1)
    nines = K.ones_like(preds) * 9
    masked = K.equal(preds - groundTruth, nines)
    return K.mean(preds - groundTruth, axis=-1)

def meanSquared(y_true, y_pred):
    groundTruth = K.argmax(y_true, axis=-1)
    preds = K.argmax(y_pred, axis=-1)
    return K.mean(K.square(preds - groundTruth), axis=-1)

def meanAbsError(y_true, y_pred):
    groundTruth = K.argmax(y_true, axis=-1)
    preds = K.argmax(y_pred, axis=-1)
    return K.mean(K.abs(preds - groundTruth), axis=-1)

cropModel = keras.models.load_model(os.path.join("models", cropModelName+".h5"), custom_objects={'meanDistance': meanDistance})

model = keras.models.load_model(os.path.join("models", modelName+".h5"), custom_objects={"meanError": meanError, "meanSquared": meanSquared, "meanAbsError": meanAbsError})

print("Assembling end-to-end model")

convNet = model.layers[1].layer
v = cropModel.output
v = Maxima2D()(v)
v = PointCrop2D(crop_size=200, mean=0.257, std=0.288, wOffset=-50, hOffset=-75)([v, cropModel.input])
v = MaxPooling2D((2,2), padding='same', name='downsampler')(v)
v = convNet(v)

visual_model = Model(inputs=cropModel.input, outputs=v)

interm_input = Input((10, numFeatures,))
layers = [l for l in model.layers]
x = layers[2](interm_input)
for i in range(3,len(layers)):
    x = layers[i](x)
linear_model = Model(inputs=interm_input, outputs=x)

logging.basicConfig(level=logging.INFO)

print("Loading session metadata")

frames = {}
softmaxs = {}
clip = VideoFileClip(os.path.join(baseDir, session, "runWisk.mp4"))
print("Duration of video (s): ",clip.duration,", FPS: ",clip.fps,", Dimensions: ",clip.size)

trackedFeatures = np.loadtxt(os.path.join(baseDir, session, "trackedFeaturesRaw.csv"), delimiter=",",skiprows=1)
mat = scipy.io.loadmat(os.path.join(baseDir, session, "runAnalyzed.mat"))

obsOnTimes = np.squeeze(mat['obsOnTimes'])
obsOnFrames = []
clipLen = len(mat['frameTimeStampsWisk'])
for i in range(len(mat['frameTimeStampsWisk'])):
    if np.isnan(mat['frameTimeStampsWisk'][i][0]):
        mat['frameTimeStampsWisk'][i][0]=-1
for i in range(len(mat['frameTimeStamps'])):
    if np.isnan(mat['frameTimeStamps'][i][0]):
        mat['frameTimeStamps'][i][0] = -1
otherFirstTime = mat['frameTimeStamps'][bisect.bisect(np.squeeze(mat['frameTimeStamps']), -1)][0]
firstTime = mat['frameTimeStampsWisk'][bisect.bisect(np.squeeze(mat['frameTimeStampsWisk']), -1)][0]
for i, obsOnTime in enumerate(obsOnTimes):
    if not (obsOnTime>min(mat['frameTimeStampsWisk'][-1], mat['frameTimeStamps'][-1]) or obsOnTime<max(firstTime, otherFirstTime)):
        obsOnFrames.append((i, bisect.bisect_left(np.squeeze(mat['frameTimeStampsWisk']), obsOnTime)))
    else:
        obsOnFrames.append((i, -1))

obsOffTimes = np.squeeze(mat['obsOffTimes'])
obsOffFrames = []
for i, obsOffTime in enumerate(obsOffTimes):
    if not (obsOffTime>min(mat['frameTimeStampsWisk'][-1], mat['frameTimeStamps'][-1]) or obsOffTime<max(firstTime, otherFirstTime)):
        obsOffFrames.append((i, bisect.bisect_left(np.squeeze(mat['frameTimeStampsWisk']), obsOffTime)))
    else:
        obsOffFrames.append((i, -1))

if len(obsOnTimes) != len(obsOffTimes):
    raise ValueError('obsOnTimes and obsOffTimes have different lengths! Please check no trials have been skipped')

temptrialFrames = list(zip(obsOnFrames, obsOffFrames))
trialFrames = []
for temp in temptrialFrames:
    if temp[0][0] != temp[1][0]:
        raise ValueError('some trial got shifted somewhere, exiting')
    trialFrames.append((temp[0][0], temp[0][1], temp[1][1]))

logging.info("Loaded features")

size = [clip.size[1], clip.size[0]]

def convertWiskStamps(wiskframe):
    time = mat['frameTimeStampsWisk'][wiskframe][0]
    return bisect.bisect_left(np.squeeze(mat['frameTimeStamps']), time)

fieldnames = ["framenum", "confidence"]
answers = [{"framenum":-1, "confidence": 0}]*len(mat['obsOnTimes'])
print("Analyzing")

for idx, framenum, endframe in tqdm(trialFrames):
    #logging.info("Obstacle on")
    #logging.info(str(framenum))
    if framenum == -1 or endframe == -1:
        continue
    findTime = 0
    initialStart = time.time()
    start = time.time()
    features = trackedFeatures[convertWiskStamps(framenum)]
    obsPos = list(map(int, [features[22], features[23]]))
    obsConf = features[24]
    nosePos = list(map(int, [features[19], features[20]]))
    image = clip.get_frame(framenum * (1 / clip.fps))
    while sum(sum(i<10 for i in image[:, (size[1]-2):, 0]))<40:
        try:
            framenum+=1
            features = trackedFeatures[convertWiskStamps(framenum)]
            obsPos = list(map(int,[features[22],features[23]]))
            obsConf = features[24]
            nosePos = list(map(int, [features[19], features[20]]))
            image = clip.get_frame(framenum * (1 / clip.fps))
        except IndexError:
            framenum = endframe+1
            break
    frameProbs = {}
    oldFrame = framenum
    if framenum >= endframe:
        print("Could not find obstacle within session. Skipping session...")
        continue
    #logging.info("Found frame!")
    #logging.info(framenum)
    findTime = time.time()-start
    loadTime = 0
    resizeTime = 0
    bottleneckTime = 0
    predictTime = 0
    cacheTime = 0
    incrementTime = 0
    needFrames = deque()
    while (nosePos[0]-obsPos[0]<50 or obsConf!=1) and framenum < endframe:
        needFrames.append(framenum)
        framenum+=1
        features = trackedFeatures[convertWiskStamps(framenum)]
        obsPos = list(map(int, [features[22], features[23]]))
        obsConf = features[24]
        nosePos = list(map(int, [features[19], features[20]]))
    for i in range(timesteps):
        needFrames.append(framenum)
        framenum+=1
    frames.clear()
    while needFrames:
        frameBatch = []
        for i in range(bs):
            if not needFrames:
                break
            frameBatch.append(needFrames.popleft())
        actualBatch = [((1.-np.reshape((clip.get_frame(i*1./clip.fps)/255.)[:,:,0], input_size+(1,)))-0.257)/0.288 for i in frameBatch]
        actualBatch = np.asarray(actualBatch)
        predBatch = visual_model.predict(actualBatch)
        frames.update({i:j for (i, j) in zip(frameBatch, predBatch)})
    framenum = oldFrame
    features = trackedFeatures[convertWiskStamps(framenum)]
    obsPos = list(map(int, [features[22], features[23]]))
    obsConf = features[24]
    nosePos = list(map(int, [features[19], features[20]]))
    needAnal = deque()
    while (nosePos[0]-obsPos[0]<50 or obsConf!=1) and framenum < endframe:
        session = []
        if framenum<timesteps:
            framenum+=1
            features = trackedFeatures[convertWiskStamps(framenum)]
            obsPos = list(map(int, [features[22], features[23]]))
            obsConf = features[24]
            nosePos = list(map(int, [features[19], features[20]]))
            continue
        if framenum>int(clip.duration*clip.fps-50):
            framenum+=1
            features = trackedFeatures[convertWiskStamps(framenum)]
            obsPos = list(map(int, [features[22], features[23]]))
            obsConf = features[24]
            nosePos = list(map(int, [features[19], features[20]]))
            continue
        for i in range(timesteps):
            frame = framenum+i
            session.append(frame)
        needAnal.append(session)
        framenum+=1
        features = trackedFeatures[convertWiskStamps(framenum)]
        obsPos = list(map(int, [features[22], features[23]]))
        obsConf = features[24]
        nosePos = list(map(int, [features[19], features[20]]))
    softmaxs.clear()
    while needAnal:
        frameBatch = []
        for i in range(bs):
            if not needAnal:
                break
            frameBatch.append(needAnal.popleft())
        actualBatch = np.asarray([[frames[i] for i in j] for j in frameBatch])
        predBatch = linear_model.predict(actualBatch)
        softmaxs.update({i[0]:j for (i,j) in zip(frameBatch, predBatch)})
    for framenum in softmaxs:
        softmax = softmaxs[framenum]
        prediction = np.argmax(softmax)
        if prediction == timesteps:
            continue
        if framenum+prediction in frameProbs:
            frameProbs[framenum+prediction]+=probDistribution[prediction]
        else:
            frameProbs[framenum+prediction]=probDistribution[prediction]
    try:
        predictedFrame = max(frameProbs.items(), key=operator.itemgetter(1))[0]
    except ValueError:
        predictedFrame=-1
    if predictedFrame != -1:
        answers[idx] = {"framenum": predictedFrame, "confidence": frameProbs[predictedFrame]}

print("Writing data")

with open(os.path.join(baseDir, sys.argv[2], "whiskerAnalyzed.csv"), 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for answer in answers:
        writer.writerow(answer)
