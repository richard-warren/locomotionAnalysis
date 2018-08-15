from typing import List, Any

import keras
from moviepy.editor import VideoFileClip
import numpy as np
import scipy.io
import csv
import sys
import os
import tensorflow as tf
import keras.backend as K

# baseDir, session
if len(sys.argv)!=3:
    print("Usage: python analyzevideo.py baseDir session")
    print(len(sys.argv))
    exit()

os.chdir(os.path.split(sys.argv[0])[0])

modelName = "dataaug4"
baseDir = sys.argv[1]
session = sys.argv[2]
fullDir = os.path.join(baseDir, session)
processor = 'CPU'

if processor == 'CPU':
    print('...analyzing paw contacts using CPU')
    config = tf.ConfigProto(device_count={'GPU':0})
    sess = tf.Session(config=config)
    K.set_session(sess)

model = keras.models.load_model(modelName+".h5")
clip = VideoFileClip(fullDir+"/runTop.mp4")
print("Duration of video(s):",clip.duration)

fieldnames = ['framenum', 'forelimb', 'hindlimb', 'notouch']

trackedFeatures = np.loadtxt(fullDir+"/trackedFeaturesRaw.csv", delimiter=",", skiprows=1)
mat = scipy.io.loadmat(fullDir+"/runAnalyzed.mat")

print("Loaded data")

def analyze(framenum, xCoord, yCoord):
    global clip
    global model
    framenum = int(framenum)
    xCoord = int(xCoord)
    yCoord = int(yCoord)
    obsPos = [xCoord, yCoord]
    size = [clip.size[1], clip.size[0]]
    image = clip.get_frame(framenum * (1 / clip.fps))
    return model.predict([np.expand_dims(np.reshape(cropFrame(image, obsPos, size)[:, :, 0],[80, 80, 1]), axis=0)/255., np.expand_dims(np.asarray([(obsPos[0]-size[1]/2)/95.4]), axis=0)])[0]

def cropFrame(image, obsPos, size):
    x = min(max(obsPos[0], 1), size[1]-2)
    y = min(max(obsPos[1], 1), size[0]-2)
    newImg = image[max(0,y-40):min(y+40, size[0]-1), max(0,x-40):min(x+40, size[1]-1)]
    if x<40:
        newImg = np.append(np.zeros((newImg.shape[0],40-x,3)).astype(int),newImg,1)
    if y<40:
        newImg = np.append(np.zeros((40-y, newImg.shape[1], 3)).astype(int), newImg, 0)
    if (size[0]-y)<41:
        newImg = np.append(newImg, np.zeros((41-size[0]+y,newImg.shape[1],3)).astype(int), 0)
    if (size[1]-x)<41:
        newImg = np.append(newImg, np.zeros((newImg.shape[0],41-size[1]+x,3)).astype(int), 1)
    return newImg

with open(fullDir+'/pawAnalyzed.csv', 'w') as csvfile:
    answers = []
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i,pos in enumerate(mat['obsPixPositions'][0]):
        features = trackedFeatures[i]
        pawX = [features[1 + i * 3] for i in range(4) if features[3 + i * 3] > 0.9]
        obsPos = list(map(int, [features[22], features[23]]))
        obsConf = features[24]
        if obsConf != 1:
            continue
        try:
            if min(abs(i - obsPos[0]) for i in pawX) > 30:
                continue
        except ValueError:
            continue
        if abs(pos-obsPos[0])<50:
            features = list(map(int,trackedFeatures[i]))
            softmax = analyze(i, obsPos[0], obsPos[1])
            answers.append({'framenum':i, 'forelimb':softmax[0], 'hindlimb':softmax[1], 'notouch':softmax[2]})
    for answer in answers:
        writer.writerow(answer)
