"""
Convert between Patient Coordinate System (PCS; Sag/Cor/Tra), Device Coordinate System (XYZ) and Gradient Coordinate System (GCS or PRS; Phase,Readout,Slice).
The coordinate systems have a common origin.

Example:
- x is a vector given in PRS-coordinates.
- prs_to_pcs() @ x = x in PCS coordinates.

Based on fantastic work from Christian Mirkes and Ali Aghaeifar.
"""

import numpy as np


internal_os = 2
pcs_directions = ['dSag', 'dCor', 'dTra']
pcs_transformations = {
    'HFS':
        [[  1,  0,  0 ],
         [  0, -1,  0 ],
         [  0,  0, -1 ]],
    'HFP':
        [[ -1,  0,  0 ],
         [  0,  1,  0 ],
         [  0,  0, -1 ]],
}


def prs_to_xyz(g):
    return pcs_to_xyz(g) @ prs_to_pcs(g)


def pcs_to_xyz(g):
    if g['patient_position'] in pcs_transformations:
        return np.array(pcs_transformations[g['patient_position']])
    else:
        raise RuntimeError("Unknown patient position")


def rps_to_xyz(g):
    return pcs_to_xyz(g) @ prs_to_pcs(g) @ rps_to_prs()


def rps_to_prs():
    return np.array(\
        [[  0,  1,  0 ],
         [  1,  0,  0 ],
         [  0,  0, -1 ]])


def prs_to_pcs(geometry):
    mat = np.eye(3)
    mat[1,1] = -1

    maindir = np.argmax(np.abs(geometry['normal']))
    if 0 == maindir:
        mat = [ [ 0, 0, 1],
                [ np.cos(geometry['angle']),  np.sin(geometry['angle']), 0],
                [-np.sin(geometry['angle']),  np.cos(geometry['angle']), 0]
              ] @ mat
        init_normal = [ 1, 0, 0 ]
    if 1 == maindir:
        mat = [ [ np.cos(geometry['angle']),  np.sin(geometry['angle']), 0],
                [ 0, 0, 1],
                [ np.sin(geometry['angle']), -np.cos(geometry['angle']), 0],
              ] @ mat
        init_normal = [ 0, 1, 0 ]
    if 2 == maindir:
        mat = [ [ np.sin(geometry['angle']), -np.cos(geometry['angle']), 0],
                [ np.cos(geometry['angle']),  np.sin(geometry['angle']), 0],
                [ 0, 0, 1],
              ] @ mat
        init_normal = [ 0, 0, 1 ]

    norm = np.linalg.norm(geometry['normal'])
    if not abs(1 - norm) < 0.001:
        raise RuntimeError(f"Normal vector is not normal?! |.| = {norm}")

    v = np.cross(init_normal, geometry['normal'])
    s = np.linalg.norm(v) # sine
    c = np.dot(init_normal, geometry['normal']) # cosine

    if s <= 0.00001:
        mat = np.eye(3) * c @ mat
    else:
        V = np.array(\
            [[     0, -v[2],  v[1] ],
             [  v[2],     0, -v[0] ],
             [ -v[1],  v[0],     0 ]])
        mat = (np.eye(3)  + V + V*V*(1-c)/s**2) @ mat

    return mat


def get_geometry(twix):
    """ Extract information about slice position from twix object
    """
    geometry = {}

    if twix['hdr']['MeasYaps']['sKSpace']['ucDimension'] == 2:
        geometry['dims'] = 2
    elif twix['hdr']['MeasYaps']['sKSpace']['ucDimension'] == 4:
        geometry['dims'] = 3
    else:
        geometry['dims'] = None

    if len(twix['hdr']['MeasYaps']['sSliceArray']['asSlice']) > 1:
        print("WARNING more than one slice. Taking first one..")

    geometry['fov'] = [
        twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]['dReadoutFOV'] * internal_os,
        twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]['dPhaseFOV'],
        twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]['dThickness']
    ]

    geometry['resolution'] = [
        twix['hdr']['MeasYaps']['sKSpace']['lBaseResolution'] * internal_os,
        twix['hdr']['MeasYaps']['sKSpace']['lPhaseEncodingLines'],
        twix['hdr']['MeasYaps']['sKSpace']['lPartitions'] if geometry['dims'] == 3 else 1
    ]

    geometry['voxelsize'] = list(np.array(geometry['fov']) / np.array(geometry['resolution']))

    geometry['normal'] = [0,0,0]
    if 'sNormal' in twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]:
        for i,d in enumerate(pcs_directions):
            geometry['normal'][i] = twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]['sNormal'].get(d, geometry['normal'][i])

    geometry['offset'] = [0,0,0]
    if 'sPosition' in twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]:
        for i,d in enumerate(pcs_directions):
            geometry['offset'][i] = twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0]['sPosition'].get(d, geometry['offset'][i])

    geometry['angle'] = twix['hdr']['MeasYaps']['sSliceArray']['asSlice'][0].get('dInPlaneRot',0)

    geometry['patient_position'] = twix['hdr']['Meas'].get('sPatPosition')

    geometry['matrix'] = rps_to_xyz(geometry).tolist()

    return geometry
