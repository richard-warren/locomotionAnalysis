import sys
import os
sys.path.append(os.path.join(os.environ['GITDIR'], 'analysis-tools', 'Python3'))
import OpenEphys

# get input arguments
file_dir = sys.argv[1]
source = sys.argv[2]
fs = sys.argv[3]
highpass = sys.argv[4]
dref = sys.argv[5]
connected_channels = sys.argv[6]

# format input arguments
fs = int(fs)
highpass = int(highpass)
if dref=='none': dref=None
connected_channels = [i=='1' for i in connected_channels]  # convert to binary vector

# pack!
print('running pack_2...')
OpenEphys.pack_2(file_dir, source=source, fs=fs, highpass=highpass, dref=dref, connected_channels=connected_channels)