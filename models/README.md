# todo
- [ ] prepare predictors for single session
- [ ] incorporate getKinematicData into autoAnalyze() and make sure diagnostic plots are produced
- [ ] single cell plots
- [ ] aggregate plots
  - [ ] sort both by peak autocorrelation AND mutual information to see if there are non-linear relationships here
  - [ ] mutual information for each cell and predictors, or cross correlations? only include high info cells in aggregate plots? does it make sense to use mutual information when model is linear? e.g. mutual info would be very high for phase predictor, but phase would be useless in model
  - [ ] use these plots to determine model transformations
- [ ] write matlab code to handle to experiment conditions (surprise, omission)

# analysis flow
- prepPredictors
  - untransformed predictors containing
    - continuous: continuous signals
    - logical: nX2 matrix of on/off times
    - event: nX1 matrix of times
- prepDesignMatrix
  - given predictors from prepPredictors, creates design matrix by
    - including only user defined predictors
    - applying transformation (e.g. convolve, raise to powers...)
- saveCellTuning
  - for each cell, compute tuning for each predictors in prepPredictors (not in design matrix)
  - for each cell-predictor combo, compute (trial X time) matrix of firing rate responses, along with x axis values, and *maybe* condition of each row (e.g. is light on, is reward omission)
    - event: true PSTHs
    - logical: interpolated 'epoch' responses
    - continuous: density estimates, along with x axis probability?
  - each predictor will need associated x axis range
  - will predictors need 'condition' labels, if i want to break down by e.g. isLightOn later?

# predictors
- continuous
  - [ ] wheel velocity
  - [X] paws (lh lf rf rh) (x y z) (position velocity)
  - [X] body angle
  - [X] whisker angle
  - [X] butt height
  - [X] jaw [along first PC]
  - [X] ear [along first PC]
  - [X] nose [along first PC]
  - [ ] satiation
  - potential additions:
    - distance to obstacle [ramping signal]
    - distance to reward [ramping signal]
- logical
  - swing stance (lh lf rf rh)
  - paw contacts (lh lf rf rh) (dorsal ventral)
  - eye open
  - obstacle on (light nolight)
  - (maybe add later) one hot vector encoding whether mouse is: running, licking, grooming, doing nothing
- event
  - whisker contact
  - licks
  - obstacle on
  - obstacle off
  - rewards

# questions
- how to handle temporal discontinuities in old sessions?
- some logical vars could be treated as events, e.g. paw contacts might better be treated as moments of contact rather than periods of contact
- how to construct distance to obstacle and distance to reward ramping signals? // want something that want signal that goes from SAME small number to 0 at whisker contact, and then fall down again // perhaps should look at all ramping cells to figure out best shape for this predictor
- how to handle moments when predictors are not known, e.g. occlusion
- is there a quantitative index that tells the non-linearity of a relationship? eg fraction of variance explained non-linearly after regressing away linear relationship?
- should maybe use continuous signal for licks instead of times, and have a 'home' position for the tongue in the mouth...

# long term todo
- [ ] figure out how to filter out poorly behaving sessions (e.g. based on velocity)
- [ ] make isSated variable, or some predictor that encodes how sated they are?
- [ ] sliding window mutual information to find optimal leads/lags for predictors, sort of like a non-linear version of cross-correlations
- [ ] should obstacle height be included somehow?
- [ ] add isBlinking based on eye tracking confidence to analyzeSession?
- [ ] housekeeping
  - [ ] make sure getKinematicData works with new analysis... will we be using this in the new project at all?
  - [ ] get rid of redundant video files
- [ ] documentation for creating training data
- [ ] add lick amplitude?
- [ ] bayesian methods for filtering tongue and whisker locations, incorporating prior information about location (in mouth, and maximally retracted)

# todo(ne)
- [X] write autoAnalyze
- [X] add omissionTimes and surpriseTimes to analyzeSession
- [X] make scoreThresh dlc vs. dpk dependent
  - [X] save metadata files for all sessions
  - [X] read files whenever confidence is needed
    - [X] fixTracking
    - [X] paw contact in analyzeSession
    - [X] getKinData
    - [X] showTracking
    - [X] showLeadingLagging
    - [X] showSingleFrameTracking
    - [X] paw contact (eddie - using low thresh for all sessions)
    - [X] wisk contact (eddie - using low thresh for all sessions)
    - [X] reanalyze all session paw contacts (without network analysis)
  - [X] add to documentation
- [X] make alignment frames automatic for problem sessions
- [X] fine tune tracking
  - [X] add to whisker training set and retrain
  - [X] add to run training set and retrain (2 lick, 8 run errors per session)
  - [X] re-analyze  and check (180917_002: 9330 15290 37888) and (200130_000)
  - [X] re-analyze all vids
  - [X] whisker contact
  - [X] body angle
  - [X] paw tracking
  - [X] paw contacts
  - [X] lick times
  - [X] check grooming on old cropped sessions AND new sessions
  - [X] whisker angle
- [X] show lick times even for low confidence frames
- [X] fix problem sessions
