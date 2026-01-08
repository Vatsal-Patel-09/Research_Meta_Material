"""
Helper utility functions.
"""

import numpy as np
from typing import Union, List


def validate_array(data: Union[List, np.ndarray]) -> np.ndarray:
    """
    Validate and convert input data to numpy array.

    Parameters
    ----------
    data : list or np.ndarray
        Input data to validate.

    Returns
    -------
    np.ndarray
        Validated numpy array.

    Raises
    ------
    ValueError
        If data is empty or contains invalid values.
    """
    arr = np.asarray(data)
    if arr.size == 0:
        raise ValueError("Input data cannot be empty")
    return arr


def normalize(data: np.ndarray) -> np.ndarray:
    """
    Normalize data to range [0, 1].

    Parameters
    ----------
    data : np.ndarray
        Input array to normalize.

    Returns
    -------
    np.ndarray
        Normalized array.
    """
    min_val = np.min(data)
    max_val = np.max(data)
    if max_val - min_val == 0:
        return np.zeros_like(data)
    return (data - min_val) / (max_val - min_val)
