import sys
sys.path.insert(0, 'D:\\github\\analysis-tools\\Python3')
from OpenEphys import pack_2



# get input arguments
file_dir = sys.argv[1]
source = sys.argv[2]
fs = sys.argv[3]
highpass = sys.argv[4]
dref = sys.argv[5]
connected_channels = sys.argv[6]

# format input arguments appropriately
fs = int(fs)
highpass = int(highpass)
if dref=='None': dref=None
connected_channels = [i=='1' for i in connected_channels] # convert to binary vector

print('running pack_2...')
pack_2(file_dir, source=source, fs=fs, highpass=highpass, dref=dref, connected_channels=connected_channels)