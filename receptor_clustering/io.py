"""CZI loading."""
import re
from pathlib import Path

import czifile
import numpy as np


def load_czi(path, channel=0, nuclear_channel=None):
    """Return (receptor stack[Z, Y, X] float32, nuclear stack or None, pixel size in um)."""
    czi = czifile.CziFile(str(path))
    arr = czi.asarray().squeeze()  # (C, Z, Y, X)
    if arr.ndim == 3:  # single channel file
        arr = arr[None]
    m = re.search(r'<Distance Id="X">\s*<Value>([^<]+)</Value>', czi.metadata())
    pixel_um = float(m.group(1)) * 1e6 if m else 1.0
    nuc = arr[nuclear_channel].astype(np.float32) if nuclear_channel is not None else None
    return arr[channel].astype(np.float32), nuc, pixel_um


def list_czi(data_dir):
    return sorted(Path(data_dir).glob("*.czi"))
