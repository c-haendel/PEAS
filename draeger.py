import os
import struct
import numpy as np
#import warnings

def load_bin(filepath, fs_target=None, max_channels_medibus=3, clamp_max_value=1e5, header_only=False):
    def mask_zero_variance(data):
        var = np.var(data, axis=0)
        zero_var_mask = var == 0
        masked = data.astype(float, copy=True)
        masked[:, zero_var_mask] = np.nan
        return masked
    try:
        with open(str(filepath), "rb") as f:
            measurements = os.path.getsize(str(filepath)) // 4358
            timestamp = np.zeros([measurements], np.float64)
            dummy = np.zeros([measurements], np.float32)
            pixel = np.zeros([measurements, 32, 32], np.float32)
            minmax = np.zeros([measurements], np.int8)
            event = np.zeros([measurements], np.int8)
            eventtxt = ["" for x in range(measurements)]
            timing = np.zeros([measurements], np.int8)
            medibus = np.zeros([measurements, 52], np.float32)

            for i in range(measurements):
                #update_progress(percent=100*i//measurements)
                timestamp[i] = 24*60*60*struct.unpack('d', f.read(8))[0] # time stamp (double) in seconds as per manual
                if(header_only is True and i == 1):
                    break
                dummy[i] = struct.unpack('f', f.read(4))[0] # dummy (float)
                for x in range(32):
                    for y in range(32):
                        pixel[i, x, y] = struct.unpack('f', f.read(4))[0] # 32*32 pixel values (float)
                minmax[i] = struct.unpack('i', f.read(4))[0] # MinMax-Flag (int)
                event[i] = struct.unpack('i', f.read(4))[0] # event marker (int)
                for k in range(30):
                    eventtxt[i] += str(struct.unpack('s', f.read(1))[0]) # 30 event text (char)
                timing[i] = struct.unpack('i', f.read(4))[0] # timing error (int)
                for k in range(52):
                    medibus[i, k] = struct.unpack('f', f.read(4))[0] # 52 medibus values (float)
                    # medibus data:
                    # 0 pressure
                    # 1 flow
                    # 2 volume ...
        if measurements <= 1:
            raise RuntimeError("File "+str(filepath)+" insufficient data points.")
        fs_bin = 1. * (measurements-1) / (timestamp[-1]-timestamp[0])

        if(header_only is False):
            if(clamp_max_value > 0):
                clamp_medibus(medibus, max_channels_medibus, clamp_max_value)

            if(not(np.all(np.isfinite(pixel)))):
                raise RuntimeError("File "+str(filepath)+" pixel is corrupt.")
            elif(not(np.all(np.isfinite(medibus)))):
                #raise RuntimeError("File "+str(filepath)+" medibus is corrupt.")
                #warnings.warn("File "+str(filepath)+" medibus is corrupt.")
                pass
            if fs_target is None: every_nth_frame=1
            else:
                every_nth_frame = int(fs_bin/fs_target)
            #if(fs_bin%fs_target>0): print("bin sampling rate", fs_bin, "not divisible by target sampling rate", fs_target)
            if(every_nth_frame > 1):
                idxs = range(0,len(pixel),every_nth_frame)
                pixel=pixel[idxs]
                medibus=medibus[idxs]
                fs_bin = fs_target
        return mask_zero_variance(pixel.astype(np.float32)), medibus.astype(np.float32), timestamp[0], fs_bin
    except(FileNotFoundError):
        raise RuntimeError("File "+str(filepath)+" not found.")

def clamp_medibus(medibus, max_channels_medibus, clamp_max_value=1e5):
    "replace large medibus entries by channel-dependent median"
    for c in range(max_channels_medibus):
        median=np.median(medibus[np.where(np.abs(medibus[:,c])<clamp_max_value)[0],c])
        medibus[np.where(np.abs(medibus[:,c])>=clamp_max_value)[0],c]=median
