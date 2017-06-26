# A direct import of matplotlib will fail if we aren't attached to a X display.
# If we can't import using the default gtk toolkit, switch to "agg" which
# doesn't require a display.

# Use agg if we know we don't have a display.
import os, sys
if('matplotlib.backends' not in sys.modules and
   (os.getenv("DISPLAY") is None or
    os.getenv("DISPLAY") == "")):
    import matplotlib
    matplotlib.use('agg')

# We still might fail, if for example we have a DISPLAY environment variable
# set for a display that isn't there anymore. So go ahead and try to import,
# but if this fails then switch to agg.
try:
    import warnings
    warnings.filterwarnings(action="ignore",
                            message=".*gdk_cursor_new_for_display: assertion `GDK_IS_DISPLAY \\(display\\)' failed.*")
    warnings.filterwarnings(action="ignore",
                            message=".*could not open display.*")
    import matplotlib.pyplot as plt
except (RuntimeError, ImportError):
    import matplotlib
    matplotlib.use('agg', warn=False, force=True)
    import matplotlib.pyplot as plt
finally:
    # Remove our filters
    warnings.filters.pop(0)
    warnings.filters.pop(0)
