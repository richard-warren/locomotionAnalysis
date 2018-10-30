import glob
import pickle
import csv


# settings
folder = 'D:\\dlcBenchmarking\\results-ricksystem\\padtests\\'
fieldnames_original = ['frame_dimensions', 'run_duration', 'nframes']
fieldnames = ['height', 'width', 'run_duration', 'nframes']




files = glob.glob(folder +  '*.pickle')



#frame_dimensions
#run_duration
#nframes


with open(folder + 'data.csv', 'w', newline='') as csvfile:
    
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    
    for i in files:
        
        infile = open(i, 'rb')
        data = pickle.load(infile)
        infile.close()
        
        run_data = {k:data['data'][k] for k in iter(fieldnames_original)}
        run_data['height'] = run_data['frame_dimensions'][0]
        run_data['width'] = run_data['frame_dimensions'][1]
        del run_data['frame_dimensions']
        
        writer.writerow(run_data)





