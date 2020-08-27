# neural recording processing pipeline
To get from a recording to a `neuralData.mat` file with instantaneous firing rates for each good unit, you must:
- spike sort (using kilosort1 and Phy)
- pick best units (using `plotQualityMetrics()`)
- prepare data for MATLAB (using `formatEphysData()`)

## spike sort
1. run `packContFiles()`
  - OpenEphys recordings have one `*.cont` file per channel, but kilosort needs all the channels 'packed' into a single `*.dat` file. `packContFiles.m` does this. It is a MATLAB wrapper for the [`pack_2`](https://github.com/richard-warren/analysis-tools/blob/c64fad3a3e11eb7221a0de6b804d864ee187ffae/Python3/OpenEphys.py#L407) function, which is part of the [`analysis-tools`](https://github.com/open-ephys/analysis-tools) library.
  - **note:** [Rick's version](https://github.com/richard-warren/analysis-tools/blob/c64fad3a3e11eb7221a0de6b804d864ee187ffae/Python3/OpenEphys.py#L407) of `pack_2` has additional functionality that we need (command average *median* filtering, and optional high pass filtering).
  - We use *no high pass filtering* at this stage of the analysis, but we do use *common median referencing*. These are the default settings for `packContFiles()`, but you can adjust them with optional Name-value arguments.


2. run kilosort1
  - Copy and paste a master_file.m from a previous ephys folder
  - Modify the folder address and parameters in both the master_file.m and corresponding config file (in D:\DCN_project\Ephys)
  - Run master_file.m located in the session ephys folder

3. run Phy (to manually curate spike kilosort1 results)
  - Run Anaconda Terminal
  - ```bash
    # navigate to the correct folder
    cd /d Z:/obstacleData/sessions/...

    # then
    phy template-gui params.py  # to initiate the Phy GUI
    ```
  - In Phy, filter KSLabel == 'good'***(only applicable in KS2)***, then sort the unit order by firing rate. for each unit, go thru the similarity view from top to bottom (ordered by similarity), do necessary merges (press g), then mark the good units as 'good' (press alt+g).
  - Criteria for merging units:
    - Two units must be on the same or adjacent channel.
    - The correlogram must not be violated.
    - The average waveform and spike shape should look very similar.

## pick best units
We found 'good' units with kilosort1, but in this final stage we pick only the *best* of those 'good' units. Also, for each good unit, we select the time window in which the unit is good. Many units have poor SNR either at the beginning or the end. This allows us to ignore those periods.

1. Run `plotQualityMetrics()`
  - This script generates quality metric plot for every good unit labelled in Phy. They are stored in `obstacleData\figures\ephys\qualityMetrics\`
  - This script doesn't need any behavior data to run.
2. Run `makeCellDataCsv()`
  - This script generate cellData.csv. Each row is a cell, and columns are 'unit_id', 'include', 'timeStart', 'timeEnd', 'location', 'notes'
  - set `include=1` for cells that meet quality standards. **todo:** what are these standards?
  - set `timeStart` and `timeEnd` to be the beginning and end of the period where the unit is usable. The units are **minutes** *with respect to the OpenEphys clock*. Note that the x axis on `plotQualityMetrics()` is also with respect to the OpenEphys clock.

## prepare for MATLAB
- Run `analyzeSession()`
  - This should be done before running the `formatEphysData`
  - This complete all neural network and other analyses for the behavioral data
- Run `formatEphysData()`
  - This creates `neuralData.mat`, which contains instantaneous firing rates and spike times (*with respect to Spike's clock*) for all of the best cells.
