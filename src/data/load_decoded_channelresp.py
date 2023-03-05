
import sys
import pickle
from scipy import io
import os
import numpy as np
import pandas as pd 
import copy


#%% Subject file
DataFile = '/Volumes/ROOT/CSNL_temp/HG/Analysis/Decoding/IEM/SJK2021/VC/channel_z_loro_min2/'
ro_path   = DataFile
ro_list   = os.listdir(ro_path)
ro_list   = sorted([s[13:17] for s in ro_list if 'sub' in s])

#%% Basic settings
with open(DataFile + "decoding_sub-" + ro_list[0] + ".pickle", 'rb') as f:
    data = pickle.load(f)