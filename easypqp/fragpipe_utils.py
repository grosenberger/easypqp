import os

import subprocess
import itertools
import numpy as np
from struct import unpack
from numpy import float32
def parse_mzml_cmd(cmd, spectralfile):
    proc = subprocess.Popen(cmd.split() + [spectralfile, str(min(os.cpu_count(), 36))], stdout=subprocess.PIPE)
    d = {}
    scannum_len = np.empty(2, dtype=np.int32)
    with proc.stdout as f:
        for _ in itertools.repeat(None, unpack('<i',f.read(4))[0]):
            f.readinto(scannum_len)
            mz_int = np.empty((2,scannum_len[1]), dtype=float32)
            f.readinto(mz_int)
            d[scannum_len[0]] = mz_int
    proc.wait()
    return d
