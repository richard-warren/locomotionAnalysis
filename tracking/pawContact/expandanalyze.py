from moviepy.editor import VideoFileClip
import numpy as np

from fastai.imports import *
from fastai.transforms import *
from fastai.conv_learner import *
from fastai.model import *
from fastai.dataset import *

import os
import shutil

import matplotlib.pyplot as plt

import csv
import sys
import scipy.io

from tqdm import tqdm

if len(sys.argv)!=3:
    print("Usage: python expandanalyze.py baseDir session")
    exit()

os.chdir(os.path.split(sys.argv[0])[0])

modelName = "168_weighted_n"
cropSize = 168
arch = resnet34
baseDir = sys.argv[1]
session = sys.argv[2]

fieldnames = ['framenum', 'fore_dorsal', 'fore_ventral', 'hind_dorsal', 'hind_ventral_high', 'hind_ventral_low', 'no_touch']

def cropFrame(image, obsPos, size, cropSize):
    cropDia = int(cropSize/2)
    x = min(max(obsPos[0], 1), size[1]-2)
    y = min(max(obsPos[1], 1), size[0]-2)
    newImg = image[max(0,y-cropDia):min(y+cropDia, size[0]-1), max(0,x-cropDia):min(x+cropDia, size[1]-1)]
    if x<cropDia:
        newImg = np.append(np.zeros((newImg.shape[0],cropDia-x,3)).astype(int),newImg,1)
    if y<cropDia:
        newImg = np.append(np.zeros((cropDia-y, newImg.shape[1], 3)).astype(int), newImg, 0)
    if (size[0]-y)<cropDia+1:
        newImg = np.append(newImg, np.zeros((cropDia+1-size[0]+y,newImg.shape[1],3)).astype(int), 0)
    if (size[1]-x)<cropDia+1:
        newImg = np.append(newImg, np.zeros((newImg.shape[0],cropDia+1-size[1]+x,3)).astype(int), 1)
    return newImg

sz = cropSize
tfms = tfms_from_model(arch, sz, aug_tfms=transforms_basic, max_zoom=1.1)
data = ImageClassifierData.from_paths("./tempgarbage", bs=1, tfms=tfms)
learn = ConvLearner.pretrained(arch, data)
learn.load(modelName)

clip = VideoFileClip(os.path.join(baseDir, session, 'run.mp4'))
print("Duration of video (s): ",clip.duration,", FPS: ",clip.fps,", Dimensions: ",clip.size)

trackedFeatures = np.loadtxt(os.path.join(baseDir, session, 'trackedFeaturesRaw.csv'), delimiter=",",skiprows=1)
mat = scipy.io.loadmat(os.path.join(baseDir, session, "runAnalyzed.mat"))

print("Loaded features")

with open(os.path.join(baseDir, session, 'pawAnalyzed.csv'), 'w') as csvfile:
    size = [168, 396]

    if not os.path.exists("./tempgarbage/"+session):
        os.makedirs("./tempgarbage/"+session)

    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    print("Preparing Data")
    for framenum, pos in enumerate(tqdm(mat['obsPixPositions'][0])):
        features = trackedFeatures[framenum]
        pawX = [features[1 + i * 3] for i in range(4) if features[3 + i * 3] > 0.9]
        obsPos = list(map(int,[features[22],features[23]]))  # todo: should not assume the obstacle position is at 22, 23
        obsConf = features[24]
        image = clip.get_frame(framenum * (1 / clip.fps))
        if obsConf != 1:
            continue
        try:
            if min(abs(i - obsPos[0]) for i in pawX) > 30:
                continue
        except ValueError:
            continue
        if not os.path.exists("./tempgarbage/"+session+"/"+str(framenum)+".png"):
            plt.imsave("./tempgarbage/"+session+"/"+str(framenum)+".png", cropFrame(image, obsPos, size, cropSize).astype(np.uint8))

    print("Running analysis")
    data = ImageClassifierData.from_paths("./tempgarbage", bs=64, tfms=tfms, test_name=session)
    learn.set_data(data)
    log_preds, y = learn.TTA(is_test=True)
    probs = np.mean(np.exp(log_preds), 0)
    answers = [{'framenum': int(os.path.split(data.test_ds.fnames[i])[1].split(".")[0]), 'fore_dorsal': j[0], 'fore_ventral': j[1], 'hind_dorsal': j[2], 'hind_ventral_high': j[3], 'hind_ventral_low': j[4], 'no_touch': j[5]} for i,j in enumerate(probs)]
    answers.sort(key=lambda x: x['framenum'])
    print("Saving data")
    for answer in answers:
        writer.writerow(answer)
    shutil.rmtree("./tempgarbage/"+session)
