# how does the brain coordinate complex behavior?

<p align="center"><img src="other/rig.png" width="500"></p>

The neural control of some simple types of movement is well understood, but how the brain coordinates complex, whole-body behavior remains a mystery. For my PhD research I developed a closed loop system in which mice run on top of a wheel and skillfully leap over motorized hurdles while I record from their brains. I developed a custom linear motion system along with custom microcontroller software that moves hurdles towards mice at the same speed that mice are running, simulating what it is like to jump over stationary objects. Using [an open source motion tracking system](https://hackaday.io/project/160744-kinemouse-wheel) I developed, we can then relate the 3D movements of mice to the activity of neurons that control these movements.

# how to...

## analyze a single session
Analyzing a session requires running four neural networks (run tracking, run contact, whisker tracking, whisker contact), and several post-processing steps to get data ready for analysis. `analyzeSession(session)` does all of this. Note that:
- you can re-analyze a session *without* re-running certain neural networks by setting any of the following optional arguments to false: `rerunRunNetwork, rerunFaceNetwork, rerunWiskContactNetwork, rerunPawContactNetwork`.
- deepposekit can analyze a run OR whisker video using the `dpkAnalysis()` wrapper
- `dpkanalysis()` can also run the implementation of deeplabcut within the deepposekit framework by specifying a deeplabcut model
- the old version of the deeplabcut analysis (for run tracking only) can be run using the `dlcAnalysis(session) `wrapper

# what is in...

## a session folder
Immediately after recording each session will contain:

file | description
--- | ---
session.smr + session.s2rx | Spike2 recording
run.mat | matlab export of Spike2 data
run.csv | metadata for run camera
wisk.csv | metadata for whisker camera
webCam.csv | metadata for webcam
run.mp4 | run video (top and bottom views concatenated)
runWisk.mp4 | whisker video
webCam.avi | webcam video
trialInfo | trial metadata (e.g. obstacle height) exported from Bonsai

`analyzeSession(session)` will create the following additional files:

file | description
--- | ---
trackedFeaturesRaw.csv | tracking for run camera
trackedFeaturesRaw_wisk.csv | tracking for whisker camera
runAnalyzed.mat | tons of useful information, e.g. frame time stamps, contact times, obstacle times, etc.
whiskerAnalyzed.csv | whisker contact analysis results
pawAnalyzed.csv | paw contact analysis results

Sessions may also include:

file | description
--- | ---
kinData.mat | "trial-ized" struct, with each row containing many metrics for a trial
ephys_* folder | contains probe recordings
neuralData.mat | processed neural data
cellData.csv | spreadsheet containing information about sorted units

## runAnalyzed.mat
todo
