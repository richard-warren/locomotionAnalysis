# todo
- [ ] fine tune tracking
  - [X] add to whisker training set and retrain
  - [X] add to run training set and retrain (2 lick, 8 run errors per session)
  - [X] re-analyze  and check (180917_002: 9330 15290 37888) and (200130_000)
  - [X] make alignment frames automatic for problem sessions
  - [X] re-analyze all vids
  - [X] whisker contact
  - [X] lick times]
  - [ ] whisker angle
  - [ ] body angle
  - [ ] fix paw contact bug (999999_999: 10612)
  - [ ] check grooming on old cropped sessions AND new sessions
  - [ ] paw tracking (180917_002: 9330 15290 37888) [should i retrain including old sessions?]
- [ ] prepare predictors for single session
- [ ] figure out how to make aggregate plots for all predictors
- [ ] sort both by peak autocorrelation AND mutual information to see if there are non-linear relationships here
- [ ] mutual information for each cell and predictors, or cross correlations? only include high info cells in aggregate plots? does it make sense to use mutual information when model is linear? e.g. mutual info would be very high for phase predictor, but phase would be useless in model
- [ ] make isSated variable, or some predictor that encodes how sated they are?
- [X] show lick times even for low confidence frames
- [X] fix problem sessions

# predictors
- continuous
  - paws (lh lf rf rh) (x y z) (position velocity)
  - paw contact (lh lf rf rh) (dorsal ventral)
  - jaw position [along first PC]
  - body angle
  - ear position [along first PC]
  - nose position [along first PC]
  - distance to obstacle [ramping signal]
  - distance to reward
  - whisker angle
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

# long term todo
- [ ] sliding window mutual information to find optimal leads/lags for predictors, sort of like a non-linear version of cross-correlations
- [ ] should obstacle height be included somehow?
- [ ] add isBlinking based on eye tracking confidence to analyzeSession?
- [ ] housekeeping
  - [ ] figure out how to deal with different scoreThresh for new sessions // which analyses will potentially be affected by this, eg paw contacts, fixTracking...
  - [ ] make sure getKinematicData works with new analysis... will we be using this in the new project at all?
  - [ ] get rid of redundant video files
- [ ] documentation for creating training data
- [ ] add lick amplitude?
