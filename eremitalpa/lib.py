"""Generic library functions."""

from collections import Counter
from typing import Callable, Iterable
import time
import functools

import numpy as np
import arviz as az
import matplotlib as mpl
import pandas as pd
import adjustText

from IPython.display import HTML


jupyter_code_cell_toggle = HTML(
    """<script>code_show = true; function code_toggle() {
        if (code_show) { $('div.input').hide(); } else {
            $('div.input').show();
        } code_show = !code_show
    } $(document).ready(code_toggle);
</script>
<font size='2'><a href='javascript:code_toggle()'>Toggle code.</a></font>
<hr />"""
)


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


def compute_errorbars(
    trace: az.InferenceData, varname: str, hdi_prob: float = 0.95
) -> np.ndarray:
    """
    Compute HDI widths for plotting with plt.errorbar.

    Args:
        trace: E.g. the output from pymc.sample.
        varname: Variable to compute error bars for.
        hdi_prob: Width of the HDI.

    Returns:
        (2, n) array of the lower and upper error bar sizes for passing to plt.errorbar.
    """
    hdi = az.hdi(trace, var_names=varname, hdi_prob=hdi_prob)[varname].values
    mean = az.extract(trace).mean(dim="sample")[varname]
    return np.stack([mean - hdi[:, 0], hdi[:, 1] - mean])


def mean(values):
    return sum(values) / len(values)


def annotate_points(
    df: pd.DataFrame, ax: mpl.axes.Axes, n: int = -1, adjust: bool = True, **kwds
):
    """
    Label (x, y) points on a matplotlib ax.

    Args:
        df: Pandas DataFrame with 2 columns, (x, y) respectively. Index contains the
         labels.
        n: Label this many points. Default (-1) annotates all points.
        ax: Matplotlib ax
        adjust: Use adjustText.adjust_text to try to prevent label overplotting.
        **kwds: Passed to adjustText.adjustText
    """
    top = df.apply(mean, axis=1).sort_values(ascending=False).head(n).index
    labels = [ax.text(*xy, site) for site, xy in df.loc[top].iterrows()]
    if adjust:
        adjustText.adjust_text(labels, ax=ax, **kwds)


def find_runs(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Find runs of consecutive items in an array.

    Args:
        arr: An array

    Returns:
        3-tuple containing run values, starts and lengths.
    """

    # ensure array
    arr = np.asanyarray(arr)
    if arr.ndim != 1:
        raise ValueError("only 1D array supported")
    n = arr.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(arr[:-1], arr[1:], out=loc_run_start[1:])
        starts = np.nonzero(loc_run_start)[0]

        # find run values
        values = arr[loc_run_start]

        # find run lengths
        lengths = np.diff(np.append(starts, n))

        return values, starts, lengths


def split_pairs(values: Iterable, separation: float = 1.0) -> list:
    """
    If values are repeated, e.g.:

        1, 5, 5, 8

    Then 'split' them by adding and subtracting half of `separation` (default=1.0 -> 0.5)
    from each item in the pair:

        1, 4.5, 5.5, 8
    """
    values = list(values)
    counts = Counter(values)
    half_separation = separation / 2
    for value, count in counts.items():
        if count == 2:
            i = values.index(value)  # index first occurrence
            j = values.index(value, i + 1)  # index second occurrence
            values[i] -= half_separation
            values[j] += half_separation
        elif count != 1:
            raise NotImplementedError(
                f"{value} occurs {count} times, only implemented pairs"
            )
    return values


def cal_months_diff(date1: pd.Timestamp, date0: pd.Timestamp) -> int:
    """Number of calendar months between two dates (date1 - date0)"""
    return (pd.Period(date1, freq="M") - pd.Period(date0, freq="M")).n
