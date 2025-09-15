########## IMPORTS ##########
import os
import scipy.io as sio
import ast
import argparse

import numpy as np
import xarray as xr
import pandas as pd

from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
from allensdk.brain_observatory.ecephys.ecephys_session import (
    EcephysSession,
    removed_unused_stimulus_presentation_columns
)

########## PARSE ARGS AND ACCESS SESSION ##########
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
type=int, help="session ID number")
ap.add_argument("-o", "--output", required=True,
help="path to output file")

args = ap.parse_args()

sessionNumber = args.input

# set path to maifest and data
manifest_path = os.path.join("D:\AllenSDK\FunctionalConnectivity", "manifest.json")
cache = EcephysProjectCache.from_warehouse(manifest=manifest_path)

# load a specific session - best is to download before running
session_id = sessionNumber #778240327 # func connec
#session_id = 715093703
session = cache.get_session_data(session_id)
print('got session')
# get  basic session info
sessionInfo = {'genotype' : session.full_genotype,
               'invalid_times' : session.invalid_times # added June '22
               }

########## PROCESS STIMULUS TABLE ##########
# get stimulus table and process required parameters
stimtable = session.get_stimulus_table();
stimtable.replace(to_replace=pd.np.nan, value="null", inplace=True)


# here we replace "null" values with None, ready to be evaled using ast.literal_eval
stimtable["pos"].replace(to_replace="null", value="None", inplace=True)
stimtable["phase"].replace(to_replace="null", value="None", inplace=True)
stimtable["size"].replace(to_replace="null", value="None", inplace=True)
stimtable["spatial_frequency"].replace(to_replace="null", value="None", inplace=True)

# convert these stimtable columns to list so we can eval them
pos = stimtable["pos"].tolist()
phase = stimtable["phase"].tolist()
size = stimtable["size"].tolist()
spatial_frequency = stimtable["spatial_frequency"].tolist()
length = len(pos)

# eval strings to become numeric, None become empty
for i in range(length):
    pos[i] = ast.literal_eval(pos[i])
    phase[i] = ast.literal_eval(phase[i])
    size[i] = ast.literal_eval(size[i])
    spatial_frequency[i] = ast.literal_eval(spatial_frequency[i])

# re-assign the lists to the stimtable dataframe
stimtable["pos"] = pos
stimtable["phase"] = phase
stimtable["size"] = size
stimtable["spatial_frequency"] = spatial_frequency

# replace any empty or "null" values with NaNs
stimtable.fillna(value=pd.np.nan, inplace=True)
stimtable.replace(to_replace="null", value=pd.np.nan, inplace=True)


########## PROCESS UNITS AND SPIKE TIMES ##########
# get generic unit info, e.g. unit location
unitInfo = session.units
unitIndex = np.array(unitInfo.index)

# dict of spike times for each unit, make matlab compatible
spikeTimes = session.spike_times
newSpikeTimes = {}
for key in spikeTimes.keys() :
    newSpikeTimes["c_" + str(key)] = [spikeTimes[key]]


########## PROCESS RUNNING AND PUPIL INFO ##########
# running speed info + pupil info
run_mid_times = session.running_speed["start_time"] + \
    (session.running_speed["end_time"] - session.running_speed["start_time"]) / 2
run_length = np.array(session.running_speed["end_time"]-session.running_speed["start_time"])
run_startTimes = np.array(session.running_speed["start_time"])
run_endTimes = np.array(session.running_speed["end_time"])
run_speed = np.array(session.running_speed["velocity"])
run_times = np.array(run_mid_times)

runInfo = {'runMidTimes' : run_times,
           'runTimeLength': run_length,
           'runStartTimes': run_startTimes,
           'runEndTimes': run_endTimes,
           'runSpeed': run_speed
           }

pupilTable = session.get_pupil_data()
if pupilTable is None: # dealing with sessions with no eye tracking
    pupilInfo = np.nan;
else:
    pupilTime = np.array(pupilTable.index)
    pupilWidth = np.array(pupilTable["pupil_width"])
    pupilHeight = np.array(pupilTable["pupil_height"])
    pupil_x = np.array(pupilTable["pupil_center_x"])
    pupil_y = np.array(pupilTable["pupil_center_y"])

    pupilInfo = {'pupilTime' : pupilTime,
              'pupilWidth': pupilWidth,
              'pupilHeight' : pupilHeight,
              'pupil_x': pupil_x,
              'pupil_y' : pupil_y
             }


########## GENERATE OUTPUT VARS DICT AND SAVE ##########
outputVars = {'sessionInfo' : sessionInfo,
              'stim_table' : {name: col.values for name, col in stimtable.items()},
              'unitInformation' : {name: col.values for name, col in unitInfo.items()},
              'unitIndex': unitIndex,
              'spikeTimes' : newSpikeTimes,
              'runInfo' : runInfo,
              'pupilInfo' : pupilInfo
            }


sio.savemat(args.output, outputVars,  long_field_names=True)
