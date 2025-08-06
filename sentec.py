import datetime
import struct
import numpy as np

def load_zri(filepath):
    def datetime_from_mltimestamp(mltimestamp):
        return(datetime.datetime.strptime("01/01/0001", "%m/%d/%Y") + datetime.timedelta(days = mltimestamp/(24*60*60*1000)))
    try:
        with open(filepath, "rb") as f:
            f.read(1) # int8
            images_list = []
            timestamps_list = []
            failing_list = []
            valid_list = []
            #failing_pattern = [0]*32
            #compensation_list = []
            #measstate_list = []
            #reconstate_list = []
            while f.read(1):
                f.seek(-1, 1)
                timestamp = datetime_from_mltimestamp(int.from_bytes(f.read(8), byteorder="little")) # uint64
                domainID = int.from_bytes(f.read(1), byteorder="little") # uint8
                numberOfDataFields = int.from_bytes(f.read(1), byteorder="little") # uint8
                failing = 0
                for i in range(numberOfDataFields):
                    dataID = int.from_bytes(f.read(1), byteorder="little")
                    payloadSize = int.from_bytes(f.read(2), byteorder="little") # payload in byte (?)
                    if domainID == 2 and dataID == 17: # extra case with differing payloadSize
                        Imagewidth = int.from_bytes(f.read(2), byteorder="little")
                        Imageheight = int.from_bytes(f.read(2), byteorder="little")
                        NumberOfMaps = int.from_bytes(f.read(2), byteorder="little")
                        f.seek(Imagewidth*Imageheight*NumberOfMaps + NumberOfMaps, 1)
            #        elif domainID == 16 and dataID == 1: # electrode quality - 2 means bad?
            #            failing_pattern_old = failing_pattern[:]
            #            for j in range(payloadSize):
            #                failing_pattern[j] = (int.from_bytes(f.read(1), byteorder="little") == 2)
            #            failing = sum(failing_pattern)
            #            compensation = (failing_pattern != failing_pattern_old)
                    elif domainID == 16 and dataID == 5: # read actual image
                        SizeWidth = int.from_bytes(f.read(1), byteorder="little")
                        SizeHeight = int.from_bytes(f.read(1), byteorder="little")
                        ZeroRef = np.empty(SizeWidth * SizeHeight)
                        for i in range(len(ZeroRef)):
                            ZeroRef[i] = struct.unpack("f", f.read(4))[0]
                        timestamps_list.append(timestamp)
                        failing_list.append(failing)
            #            compensation_list.append(compensation)
                        if ZeroRef[16*32+16] == 1: # do not include images with negative one in center
                            images_list.append(np.full((SizeWidth, SizeHeight), np.nan))
                            valid_list.append(False)
                        else:
                            images_list.append(-ZeroRef.reshape(SizeWidth, SizeHeight).copy())
                            valid_list.append(True)
            #        elif domainID == 16 and dataID == 10: # ReconState und MeasState
            #            ReconState = int.from_bytes(f.read(1), byteorder="little")
            #            reconstate_list.append(ReconState)
            #            MeasState = int.from_bytes(f.read(1), byteorder="little")
            #            measstate_list.append(MeasState)
                    else: # discard all other info
                        #print(int.from_bytes(f.read(payloadSize), byteorder="little"))
                        f.seek(payloadSize, 1)
        if len(images_list) == 0:
            raise RuntimeError(f"File {filepath} has insufficient data points.")
        return(np.array(images_list)[:,:,:], timestamps_list, failing_list, valid_list)
    except(FileNotFoundError):
        raise RuntimeError(f"File {filepath} not found, returning empty ImageData Object.")
