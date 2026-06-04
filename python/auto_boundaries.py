"""
Automatic Boundary Detection via Diagonal Profile Analysis

Identifies bicluster boundaries in a reordered similarity matrix by:
1. Computing a diagonal profile (local neighborhood mean along the diagonal)
2. Finding valleys (local minima) using scipy peak detection

Used by BiSer to segment the TSP-reordered similarity matrix into biclusters.
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema, find_peaks


def find_auto_boundaries(df_or_matrix, valley, smooth=0, window=10,
                         sigma=3, order=5, prominence=0.01, distance=5):
    """
    Detect bicluster boundaries from a reordered similarity matrix.

    Parameters
    ----------
    df_or_matrix : numpy array or pandas DataFrame
        Square similarity matrix (reordered by TSP).
    valley : str
        Method for valley detection:
        - 'find_peaks': scipy find_peaks on inverted profile (recommended)
        - 'argrelextrema': scipy argrelextrema with elbow selection
    smooth : str or int
        If 'gaussian', apply Gaussian smoothing before peak detection.
    window : int
        Window size for diagonal profile computation.
    sigma : float
        Sigma for Gaussian smoothing (only used if smooth='gaussian').
    order : int
        Order parameter for argrelextrema (only used if valley='argrelextrema').
    prominence : float
        Minimum prominence for find_peaks (only used if valley='find_peaks').
    distance : int
        Minimum distance between peaks for find_peaks.

    Returns
    -------
    list of int
        Indices of detected bicluster boundaries.
    """
    window = int(window)
    order = int(order)

    # Convert DataFrame to numpy array if necessary
    if hasattr(df_or_matrix, "to_numpy"):
        matrix = df_or_matrix.to_numpy()
    else:
        matrix = df_or_matrix

    # Compute diagonal profile: mean of local neighborhood along the diagonal
    def diagonal_profile(matrix, window=10):
        n = matrix.shape[0]
        profile = []
        for i in range(n):
            lo = max(0, i - window)
            hi = min(n, i + window + 1)
            vals = [matrix[i, j] for j in range(lo, hi)]
            profile.append(np.mean(vals))
        return np.array(profile)

    profile = diagonal_profile(matrix, window=window)

    # Optional smoothing
    if smooth == 'gaussian':
        profile_smooth = gaussian_filter1d(profile, sigma=sigma)
    else:
        profile_smooth = profile

    # Valley detection
    if valley == 'argrelextrema':
        local_minima = argrelextrema(profile_smooth, np.less, order=order)[0]
        valley_depths = profile_smooth[local_minima]

        if len(valley_depths) < 2:
            return []

        sorted_indices = np.argsort(valley_depths)
        sorted_depths = valley_depths[sorted_indices]
        depth_diffs = np.diff(sorted_depths)

        elbow_index = np.argmax(depth_diffs) + 1
        selected_valley_indices = sorted_indices[:elbow_index]
        auto_boundaries = sorted(local_minima[selected_valley_indices])

    elif valley == 'find_peaks':
        inverted = -profile_smooth
        local_minima, properties = find_peaks(inverted,
                                               prominence=prominence,
                                               distance=distance)
        auto_boundaries = local_minima

    else:
        raise ValueError(f"Unknown valley method: {valley}")

    return auto_boundaries
