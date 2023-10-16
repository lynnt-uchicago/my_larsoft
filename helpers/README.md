# Python Modules/Helper Functions

This directory contains a collection of Python modules. While these python modules can be accessed by a Python script/Jupyter notebook from any directory with an absolute path, it can be useful to place these modules in a directory that is in your `PYTHONPATH` environment variable. This will allow you to import these modules from any directory without having to specify the absolute path to the module. If you are using a virtual environment, it can also be placed inside the `site-packages` directory of your virtual environment.

## `nue` module

This module contains a collection of helper functions for performing the nue selection on processed caf files. In order to go from a regular flat caf file (the output of running `cafmakerjob<...>.fcl` in LArSoft), to the processed caf file, you need to run Gray's workflow, details available [here](https://github.com/SBNSoftware/sbnana/tree/develop/sbnana/SBNAna/icarus-analysis-villiage/pyana/dimuon-tools), all glory to Gray. This may seem overkill for a handful of caf files, but is extremely necessary when performing a selection over an entire production.

These helper functions may be helpful for any selectio, but there some assumptions/dependencies that are specific to the nue selection.

To use this module you can import it as such:

``` python
import nue as nue
# this is useful if you are making changes to nue.py and want to reload the module without restarting the kernel (in jupyter notebook)
importlib.reload(nue)
...
# example calling of the flatten_df function from the nue module 
flat_df = nue.flatten_df(df) 
```

### Notes

- This module contains some hard-coded fixes for bugs in the cafs, such as the pfp track length and the bestPlane for a shower.

## `wc_helper` module

This module contains helper functions for performing WireCell validation and analysis. The functions are mainly used on the output of `celltree_<exp>.fcl` and/or `dump_waveform.C`.
