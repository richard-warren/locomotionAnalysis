import WhiskiWrap
#import multiprocessing
import tables
import pandas
import datetime
import os



if __name__ == '__main__':
	
	start = datetime.datetime.now()

	video_filename = 'C:/Users/rick/Google Drive/columbia/obstacleData/sessions/171231_002/runWisk.mp4'
	#video_filename = "C:/Users/rick/Anaconda3/envs/whiskiWrapEnv/Scripts/WhiskiWrap/test_video2.mp4"
	input_reader = WhiskiWrap.FFmpegReader(video_filename, duration=None, write_stderr_to_screen=False)
	tiffs_to_trace_directory = 'tiffs'
	WhiskiWrap.interleaved_reading_and_tracing(input_reader, tiffs_to_trace_directory, chunk_size=200, delete_tiffs=True, n_trace_processes=4, verbose=True, stop_after_frame=None) 

	# Load results from hdf5 file
	with tables.open_file('whiskers.h5') as fi:                               
		hdf5_whiskers = pandas.DataFrame.from_records(fi.root.summary.read())


	stop = datetime.datetime.now()
	print "Time required: %s" % (stop - start)