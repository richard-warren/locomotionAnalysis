# locomotionAnalysis

<p align="center"><img src="other/rig.png" width="500"></p>

The neural control of some simple types of movement is well understood, but how the brain coordinates complex, whole-body behavior remains a mystery. For my PhD research I developed a closed loop system in which mice run on top of a wheel and skillfully leap over motorized hurdles while I record from their brains. I developed a custom linear motion system along with custom microcontroller software that moves hurdles towards mice at the same speed that mice are running, simulating what it is like to jump over stationary objects. Using [an open source motion tracking system](https://hackaday.io/project/160744-kinemouse-wheel) I developed, we can then relate the 3D movements of mice to the activity of neurons that control these movements.

# how to

#### analyze a single session
Analyzing a session requires running four neural networks (run tracking, run contact, whisker tracking, whisker contact), and several post-processing steps to get data ready for analysis. `analyzeSession(session)` does all of this. Note that:
- you can re-analyze a session *without* re-running certain neural networks by setting any of the following optional arguments to false: `rerunRunNetwork, rerunFaceNetwork, rerunWiskContactNetwork, rerunPawContactNetwork`.
- deepposekit can analyze a run OR whisker video using the `dpkAnalysis()` wrapper
- `dpkanalysis()` can also run the implementation of deeplabcut within the deepposekit framework by specifying a deeplabcut model
- the old version of the deeplabcut analysis (for run tracking only) can be run using the `dlcAnalysis(session) `wrapper

#### what is in a session folder?
Among other things, we have:

file | description
--- | ---
run.csv | camera metadata for run camera
wisk.csv | camera metadata for wisk camera
pawAnalyzed.csv | paw contact analysis results
whiskerAnalyzed.csv | whisker contact analysis results
trackedFeaturesRaw.csv | tracking for run camera
trackedFeaturesRaw_wisk.csv | tracking for wisk camera
