# ruff: noqa :E402
from __future__ import annotations

from pathlib import Path

from spinifex import h5parm_tools
from spinifex.vis_tools import ms_tools

"""example how to get rm values for an ms and write them to an h5parm file"""
### Required to load local data for example - not needed for normal use
from importlib import resources

with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
    spinifex_data = datapath
# We set these options to use the data packaged with spinifex
# Unsetting them will cause the function to download the data from the internet
iono_kwargs = {"prefix": "cod", "output_directory": spinifex_data}
###
ms = "test.MS"
# Unpacking of local test.MS - not needed for normal use

import shutil

with resources.as_file(resources.files("spinifex.data.tests")) as test_data:
    zipped_ms = test_data / "test.ms.zip"
    shutil.unpack_archive(zipped_ms, "./")
#
# create a Path object for ms
ms_path = Path(ms)
# get metadata from ms, needed to get station_names
ms_metadata = ms_tools.get_metadata_from_ms(ms_path)
# get a dictionary with rm objects, keys are the station names
rms = ms_tools.get_rm_from_ms(
    ms_path, iono_kwargs=iono_kwargs, use_stations=ms_metadata.station_names
)

h5parm_name = "test.h5"
# write to an h5parm object, this can be a new or existing h5parm
h5parm_tools.write_rm_to_h5parm(rms=rms, h5parm_name=h5parm_name)
