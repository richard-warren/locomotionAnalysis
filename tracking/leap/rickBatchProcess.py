# Yuck -- hackish relative importing:
import os, sys
print(os.getcwd())
if os.path.exists(os.path.join(os.getcwd(), "leap", "predict_box.py")):
    leap_dir = os.path.join(os.getcwd())
elif os.path.exists(os.path.join(os.path.dirname(os.getcwd()), "leap", "predict_box.py")):
    leap_dir = os.path.join(os.path.dirname(os.getcwd()))
# leap_dir = ".." # replace this with the absolute path if imports are not working
sys.path.append(leap_dir) # add path to repository root

import numpy as np
import cv2
import h5py
from time import time

import keras
import keras.models
from leap.predict_box import convert_to_peak_outputs
from leap.utils import versions

versions(list_devices=True)







# Media file path
video_path = "C:\\Users\\rick\Desktop\\github\\leap\\leap\\batchTest.mp4"

# Trained network path
model_path = "C:\\Users\\rick\\Desktop\\github\\leap\\models\\BermanFlies\\FlyAging-DiegoCNN_v1.0_filters=64_rot=15_lrfactor=0.1_lrmindelta=1e-05_01\\final_model.h5"

# Predictions output path
save_path = "C:\\Users\\rick\\Desktop\\github\\leap\\leap\\batchTest.mp4\\072212_163153.preds.h5"

# Number of frames to read before predicting (higher = faster, but limited by RAM)
chunk_size = 5000

# Number of frames to evaluate at once on the GPU (higher = faster, but limited by GPU memory)
batch_size = 16








t0_all = time()

# Load model and convert to peak-coordinate output
model = convert_to_peak_outputs(keras.models.load_model(model_path))
print("Model:", model_path)
print("    Input:", str(model.input_shape))
print("    Output:", str(model.output_shape))

# model = keras.utils.multi_gpu_model(model, gpus=2)

# Open video for reading
reader = cv2.VideoCapture(video_path)
num_samples = int(reader.get(cv2.CAP_PROP_FRAME_COUNT))

# Initialize
positions_pred = []
conf_pred = []
buffer = []
samples_predicted = 0
reading_runtime = 0
prediction_runtime = 0
done = False

# Process video chunk-by-chunk
while not done:
    t0_reading = time()
    # Read and finish if no frame was retrieved
    returned_frame, I = reader.read()
    done = not returned_frame
    reading_runtime += time() - t0_reading
    
    # Add current frame to buffer
    if not done:
        buffer.append(I[...,0])
    
    # Do we have anything to predict?
    if len(buffer) >= chunk_size or (done and len(buffer) > 0):
        t0_prediction = time()
        
        # Predict on buffer
        Y = model.predict(np.stack(buffer, axis=0)[...,None], batch_size=batch_size)
        
        # Save
        positions_pred.append(Y[:,:2,:].astype("int32"))
        conf_pred.append(Y[:,2,:].squeeze())
        
        # Empty out buffer container
        buffer = []
        
        # Performance stats
        samples_predicted += len(Y)
        prediction_runtime += time() - t0_prediction
        elapsed = time() - t0_all
        fps = samples_predicted / elapsed
        print("Predicted: %d/%d frames | Elapsed: %.1f min / %.1f FPS / ETA: %.1f min" %
              (samples_predicted, num_samples, elapsed / 60, fps, (num_samples - samples_predicted) / fps / 60))
        
# Close video reader
reader.release()

# Merge arrays
positions_pred = np.concatenate(positions_pred, axis=0)
conf_pred = np.concatenate(conf_pred, axis=0)

# Report performance stats
print("Finished predicting %d frames." % samples_predicted)
print("    Prediction | Runtime: %.2f min / %.3f FPS" % (prediction_runtime / 60, samples_predicted / prediction_runtime))
print("    Reading    | Runtime: %.2f min / %.3f FPS" % (reading_runtime / 60, samples_predicted / reading_runtime))

# Save
if os.path.exists(save_path):
    os.remove(save_path)
with h5py.File(save_path, "w") as f:
        f.attrs["num_samples"] = num_samples
        f.attrs["video_path"] = video_path
        f.attrs["model_path"] = model_path

        ds_pos = f.create_dataset("positions_pred", data=positions_pred, compression="gzip", compression_opts=1)
        ds_pos.attrs["description"] = "coordinate of peak at each sample"
        ds_pos.attrs["dims"] = "(sample, [x, y], joint) === (sample, [column, row], joint)"

        ds_conf = f.create_dataset("conf_pred", data=conf_pred, compression="gzip", compression_opts=1)
        ds_conf.attrs["description"] = "confidence map value in [0, 1.0] at peak"
        ds_conf.attrs["dims"] = "(sample, joint)"

        total_runtime = time() - t0_all
        f.attrs["reading_runtime_secs"] = reading_runtime
        f.attrs["prediction_runtime_secs"] = prediction_runtime
        f.attrs["total_runtime_secs"] = total_runtime
        
    
print("Saved:", save_path)

print("Total runtime: %.1f mins" % (total_runtime / 60))
print("Total performance: %.3f FPS" % (samples_predicted / total_runtime))