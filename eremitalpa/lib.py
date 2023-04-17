"""Generic library functions."""

from typing import Callable
import time
import functools

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
