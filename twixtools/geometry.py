#!/usr/bin/python
"""Extract information to convert between
Patient Coordinate System (PCS; Sag/Cor/Tra), Device Coordinate System (XYZ) and Gradient Coordinate System (GCS or PRS; Phase,Readout,Slice).
Example:
- x is a vector given in PRS-coordinates.
- prs_to_pcs() @ x = x in PCS coordinates.

Based on work by Christian Mirkes and Ali Aghaeifar.
"""


import argparse
import json

import numpy as np

import twixtools


internal_os = 2
pcs_directions = ["dSag", "dCor", "dTra"]

# p. 418 - pcs to dcs
pcs_transformations = {
    "HFS": [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    "HFP": [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    "FFS": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
}


def rps_to_prs():
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])


def get_plane_orientation(geometry):
    if not abs(1 - np.linalg.norm(geometry["normal"])) < 0.001:
        raise RuntimeError(f"Normal vector is not normal: |x| = {norm}")

    maindir = np.argmax(np.abs(geometry["normal"]))
    if 0 == maindir:
        mat = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]  # @ mat // inplane mat
    if 1 == maindir:
        mat = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    if 2 == maindir:
        mat = np.eye(3)

    init_normal = np.zeros(3)
    init_normal[maindir] = 1

    v = np.cross(init_normal, geometry["normal"])
    s = np.linalg.norm(v)
    c = np.dot(init_normal, geometry["normal"])

    if s <= 0.00001:
        mat = np.eye(3) * c @ mat
    else:
        V = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        mat = (np.eye(3) + V + V * V * (1 - c) / s ** 2) @ mat

    return mat


def get_inplane_rotation(geometry):
    mat = [
        [-np.sin(geometry["angle"]), np.cos(geometry["angle"]), 0],
        [-np.cos(geometry["angle"]), -np.sin(geometry["angle"]), 0],
        [0, 0, 1],
    ]
    return np.array(mat)


def prs_to_pcs(geometry):
    mat = get_inplane_rotation(geometry)
    mat = get_plane_orientation(geometry) @ mat
    return mat


def pcs_to_xyz(g):
    if g["patient_position"] in pcs_transformations:
        return np.array(pcs_transformations[g["patient_position"]])
    else:
        raise RuntimeError(f"Unknown patient position: {g['patient_position']}")


def prs_to_xyz(g):
    return pcs_to_xyz(g) @ prs_to_pcs(g)


def rps_to_xyz(g):
    return prs_to_xyz(g) @ rps_to_prs()


def get_geometry(twix):
    geometry = {}

    if twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 2:
        geometry["dims"] = 2
    elif twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 4:
        geometry["dims"] = 3
    else:
        geometry["dims"] = None

    if len(twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"]) > 1:
        print("WARNING more than one slice. Taking first one..")

    geometry["fov"] = [
        twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dReadoutFOV"]
        * internal_os,
        twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dPhaseFOV"],
        twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dThickness"],
    ]

    geometry["resolution"] = [
        twix["hdr"]["MeasYaps"]["sKSpace"]["lBaseResolution"] * internal_os,
        twix["hdr"]["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
        twix["hdr"]["MeasYaps"]["sKSpace"]["lPartitions"]
        if geometry["dims"] == 3
        else 1,
    ]

    geometry["voxelsize"] = list(
        np.array(geometry["fov"]) / np.array(geometry["resolution"])
    )

    geometry["normal"] = [0, 0, 0]
    if "sNormal" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]:
        for i, d in enumerate(pcs_directions):
            geometry["normal"][i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][
                0
            ]["sNormal"].get(d, geometry["normal"][i])

    geometry["offset"] = [0, 0, 0]
    if "sPosition" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]:
        for i, d in enumerate(pcs_directions):
            geometry["offset"][i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][
                0
            ]["sPosition"].get(d, geometry["offset"][i])

    geometry["angle"] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0].get(
        "dInPlaneRot", 0
    )

    if "tPatientPosition" in twix["hdr"]["Meas"]:
        geometry["patient_position"] = twix["hdr"]["Meas"].get("tPatientPosition")
    elif "sPatPosition" in twix["hdr"]["Meas"]:
        geometry["patient_position"] = twix["hdr"]["Meas"].get("sPatPosition")
    else:
        geometry["patient_position"] = None

    geometry["matrix"] = rps_to_xyz(geometry).tolist()

    return geometry


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("raw", help="Raw data file")
    parser.add_argument(
        "outfile", nargs="?", default="-", help=".json file to store the output"
    )
    args = parser.parse_args()

    twix = twixtools.read_twix(args.raw, parse_data=False)
    geometry = get_geometry(twix[-1])
    if args.outfile == "-":
        print(json.dumps(geometry, indent=4, sort_keys=True))
    else:
        with open(args.outfile, "w") as f:
            json.dump(geometry, f)
