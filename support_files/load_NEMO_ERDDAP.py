# Copyright 2016 The Salish Sea NEMO Project and
# The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for Salish Sea NEMO postprocessing and visualization
"""

import numpy                as np
import xarray               as xr
import os
from mpl_toolkits.basemap import Basemap

def load_NEMO(
        timerange, depth, w_grid=[100, 398, 230, 570],
        ERDDAP_url='https://salishsea.eos.ubc.ca/erddap/griddap'
):
    """Loads and processes Salish Sea NEMO 3.4 salinity and velocity fields
    on the specified timerange, depth, domain window, and vector spacing.
    """
    
    # Slice objects for xarray indexing
    tslice = slice(timerange[0], timerange[1])
    xslice = slice(w_grid[0], w_grid[1])
    yslice = slice(w_grid[2], w_grid[3])
    
    # Load Salish Sea NEMO results from ERDDAP server as xarray DataSets
    NEMO_grid = xr.open_dataset(os.path.join(ERDDAP_url, 'ubcSSnBathymetry2V1'))
    NEMO_trc  = xr.open_dataset(os.path.join(ERDDAP_url, 'ubcSSn3DTracerFields1hV1'))
    NEMO_u    = xr.open_dataset(os.path.join(ERDDAP_url, 'ubcSSn3DuVelocity1hV1'))
    NEMO_v    = xr.open_dataset(os.path.join(ERDDAP_url, 'ubcSSn3DvVelocity1hV1'))
    
    # Extract results as xarray DataArrays
    lon      = NEMO_grid.longitude.sel(gridX=xslice, gridY=yslice)
    lat      = NEMO_grid.latitude.sel( gridX=xslice, gridY=yslice)
    salinity = NEMO_trc.salinity.sel(  gridX=xslice, gridY=yslice,
                                time=tslice).sel(depth=depth, method='nearest')
    u        = NEMO_u.uVelocity.sel(   gridX=xslice, gridY=yslice,
                                time=tslice).sel(depth=depth, method='nearest')
    v        = NEMO_v.vVelocity.sel(   gridX=xslice, gridY=yslice,
                                time=tslice).sel(depth=depth, method='nearest')
    
    # Unstagger currents
    u = np.add(u[..., :-1   ], u[..., 1:   ]) / 2
    v = np.add(v[..., :-1, :], v[..., 1:, :]) / 2
    v = v.reindex_like(u)
    
    # Rotate currents
    theta_rad = 29 * np.pi / 180
    u_map = u * np.cos(theta_rad) - v * np.sin(theta_rad)
    v_map = u * np.sin(theta_rad) + v * np.cos(theta_rad)
    
    # Return as dict
    NEMO = {'lon': lon, 'lat': lat, 'salinity': salinity, 'u': u_map, 'v': v_map}
    
    return NEMO


def plot_map(w_map, projection='lcc', resolution='h'):
    """Plot a Basemap instance on the given projection and resolution
    with map boundaries defined by w_map
    """
    
    # Make projection
    m = Basemap(projection=projection, resolution=resolution,
                lon_0=(w_map[1] - w_map[0]) / 2 + w_map[0],
                lat_0=(w_map[3] - w_map[2]) / 2 + w_map[2],
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3])

    # Add features and labels
    m.drawcoastlines()
    m.fillcontinents(color='burlywood')
    m.drawmeridians(np.arange(w_map[0], w_map[1], 0.5), labels=[0, 0, 0, 1])
    m.drawparallels(np.arange(w_map[2], w_map[3], 0.5), labels=[1, 0, 0, 0])
    
    return m
