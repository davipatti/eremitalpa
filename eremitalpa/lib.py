"""Generic library functions."""

from typing import Callable
import time
import functools


def log_df_func(f: Callable, *args, **kwargs):
    """
    Callable should return a DataFrame. Report time taken to call a function, and the
    shape of the resulting DataFrame.
    """

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        tick = time.time()
        out = f(*args, **kwargs)
        tock = time.time()
        print(f"{f.__name__} {tock - tick:.3f}s. shape: ({out.shape})")
        return out

    return wrapped
