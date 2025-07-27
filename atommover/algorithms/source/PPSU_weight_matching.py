# Authors:
# source code in C (bttlThreshold): [PPSU2023](https://inria.hal.science/hal-04146298)
# wrapper: ChatGPT

import ctypes
import os
import subprocess
import sys

# Locate the shared library in the c_code directory
BASE_DIR = os.path.abspath(os.path.join(__file__, "../../../.."))

# Path to the PPSU2023 directory (in the project root)
PPSU_DIR = os.path.join(BASE_DIR, "PPSU2023")

# Path to the shared library and setup script
LIB_PATH = os.path.join(PPSU_DIR, "libmatching_for_PPSU.so")
SETUPC_PATH = os.path.join(PPSU_DIR, "setupc.py")

def build_shared_library():
    """Run PPSU2023/setupc.py to build the shared library."""
    if not os.path.isfile(SETUPC_PATH):
        raise RuntimeError(f"setupc.py not found at expected path: {SETUPC_PATH}")
    print(f"[INFO] Attempting to build shared library via {SETUPC_PATH}...")
    subprocess.run([sys.executable, SETUPC_PATH, "build_ext"], check=True, cwd=PPSU_DIR)
    print("[INFO] Build completed.")

try:
    lib = ctypes.CDLL(LIB_PATH)
except OSError as e:
    print(f"[WARNING] Failed to load shared library at {LIB_PATH}: {e}")
    build_shared_library()
    lib = ctypes.CDLL(LIB_PATH)

libc = ctypes.CDLL(None)
stdout_fileno = sys.stdout.fileno()
libc.fflush(None)  # Flush C stdio buffers

# Define the function signature for bttlThreshold
lib.bttlThreshold.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # col_ptrs
    ctypes.POINTER(ctypes.c_int),  # col_ids
    ctypes.POINTER(ctypes.c_double),  # col_vals
    ctypes.c_int,  # n
    ctypes.c_int,  # m
    ctypes.POINTER(ctypes.c_int),  # match
    ctypes.POINTER(ctypes.c_int),  # row_match
    ctypes.POINTER(ctypes.c_int),  # row_ptrs
    ctypes.POINTER(ctypes.c_int),  # row_ids
    ctypes.POINTER(ctypes.c_double),  # row_vals
    ctypes.POINTER(ctypes.c_int),  # fend_cols
    ctypes.POINTER(ctypes.c_int),  # fend_rows
    ctypes.c_int,  # lbapAlone
    ctypes.POINTER(ctypes.c_double),  # thrshld_g
    ctypes.c_int,  # sprankknown
]
lib.bttlThreshold.restype = ctypes.c_int  # Returns the number of iterations

def bttl_threshold(col_ptrs, col_ids, col_vals, n, m, sprankknown=0, lbapAlone=1):
    """
    Python wrapper for the bttlThreshold function in the shared library.

    Args:
        col_ptrs (list[int]): Column pointers (CSR format).
        col_ids (list[int]): Row indices (CSR format).
        col_vals (list[float]): Edge weights in CSR format.
        n (int): Number of columns.
        m (int): Number of rows.
        sprankknown (int): Structural rank of the matrix (default: 0).

    Returns:
        dict: Matching results, including column-to-row, row-to-column mappings, and threshold.
    """
    # Convert inputs to ctypes
    col_ptrs = (ctypes.c_int * len(col_ptrs))(*col_ptrs)
    col_ids = (ctypes.c_int * len(col_ids))(*col_ids)
    col_vals = (ctypes.c_double * len(col_vals))(*col_vals)
    match = (ctypes.c_int * n)(-1)  # Initialize match array with -1
    row_match = (ctypes.c_int * m)(-1)  # Initialize row_match array with -1
    row_ptrs = (ctypes.c_int * (m + 1))()  # Placeholder for row pointers
    row_ids = (ctypes.c_int * len(col_ids))()  # Placeholder for row indices
    row_vals = (ctypes.c_double * len(col_vals))()  # Placeholder for row values
    fend_cols = (ctypes.c_int * n)()  # Placeholder for fend_cols
    fend_rows = (ctypes.c_int * m)()  # Placeholder for fend_rows
    thrshld_g = ctypes.c_double()  # Threshold value

    # Call the C function
    iterations = lib.bttlThreshold(
        col_ptrs, col_ids, col_vals, ctypes.c_int(n), ctypes.c_int(m),
        match, row_match, row_ptrs, row_ids, row_vals,
        fend_cols, fend_rows, ctypes.c_int(lbapAlone), ctypes.byref(thrshld_g), sprankknown
    )

    # Return results
    return {
        "iterations": iterations,
        "match": list(match),
        "row_match": list(row_match),
        "threshold": thrshld_g.value,
    }
