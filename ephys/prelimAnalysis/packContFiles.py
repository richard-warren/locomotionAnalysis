import sys
sys.path.insert(0, 'D:\\DCN_project\\Github\\analysis-tools\\Python3')
from OpenEphys import pack_2



# get input arguments
file_dir = 'Z:\\obstacleData\\sessions\\200703_000\\ephys_2020-07-03_13-48-14'
source = 107
fs = 30000
highpass = 0
dref = 'med'
connected_channels = 'all'

# format input arguments appropriately
fs = int(fs)
highpass = int(highpass)
if dref=='None': dref=None


print('running pack_2...')
pack_2(folderpath = file_dir, filename = '', channels = 'all', chprefix = 'CH',
           dref = dref, session = '0', source = '107')