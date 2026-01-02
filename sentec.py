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


def load_eit(filepath):
    dt0 = datetime.datetime.strptime("01/01/0001", "%m/%d/%Y")
    dt_from_ms = lambda ms: dt0 + datetime.timedelta(milliseconds=int(ms))

    def r(f, fmt):
        return struct.unpack(fmt, f.read(struct.calcsize(fmt)))[0]

    def rvec_i32(f, n):
        return np.frombuffer(f.read(4*n), dtype="<i4", count=n)

    try:
        with open(filepath, "rb") as f:
            fmtver = r(f, "<i")
            if fmtver not in (4, 5):
                raise RuntimeError(f"Unsupported SenTec/Swisstom .eit format_version={fmtver}")
            hdrsz = r(f, "<i")
            f.seek(16, 1)
            f.seek(hdrsz, 0)

            # Landquart LQ4/LQ5 constants (EIDORS reader)
            EOFF, IQN, VIN, POSN = 328, 2048, 64, 3
            amp = 2.048 / (2**20 * 360 * 1000)  # "simple guess" scaling as implemented in EIDORS

            iq_cols, vi_cols, pos_cols, tabs, trel, evts = [], [], [], [], [], []

            while True:
                b = f.read(EOFF)
                if len(b) < EOFF:
                    break
                f.seek(-EOFF, 1)

                tAbs_ms = r(f, "<q")
                ft = r(f, "<i")
                pl = r(f, "<i")

                if ft == 1:
                    ev = {"timestamp_ms": tAbs_ms}
                    if pl >= 4:
                        ev["eventId"] = r(f, "<i")
                        pl -= 4
                    evts.append(ev)
                    if pl > 0:
                        f.seek(pl, 1)
                    continue

                if ft != 0:
                    if pl > 0:
                        f.seek(pl, 1)
                    continue

                _hdr = rvec_i32(f, 15)
                trel.append(int(_hdr[4]))
                pos_cols.append(rvec_i32(f, POSN))
                vi_cols.append(rvec_i32(f, VIN))
                iq_cols.append(rvec_i32(f, IQN))
                tabs.append(dt_from_ms(tAbs_ms))

                skip = pl - (4*IQN + EOFF)
                if skip > 0:
                    f.seek(skip, 1)

            if not iq_cols:
                raise RuntimeError(f"File {filepath} has insufficient data points.")

            iq = np.stack(iq_cols, axis=1)
            vi = np.stack(vi_cols, axis=1)
            pos = np.stack(pos_cols, axis=1) if pos_cols else np.zeros((POSN, iq.shape[1]), dtype=np.int32)

            vv = amp * (iq[0::2].astype(np.float64) + 1j * iq[1::2].astype(np.float64))
            elecImps = vi[0::2].astype(np.float64) + 1j * vi[1::2].astype(np.float64)

            return vv, tabs, trel, elecImps, pos, evts, {"format_version": fmtver, "header_size": hdrsz}

    except FileNotFoundError:
        raise RuntimeError(f"File {filepath} not found, returning empty EITData Object.")
