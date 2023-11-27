#!/usr/bin/python
import numpy as np

internal_os = 2
pcs_directions = ["dSag", "dCor", "dTra"]

# p. 418 - pcs to dcs
pcs_transformations = {
    "HFS": [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    "HFP": [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    "FFS": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
}


class Geometry:
    """Get geometric information from twix dict

    During initialization, information about slice geometry is copied from the supplied twix dict.
    Methods for conversion between the different coordinate systems
    Patient Coordinate System (PCS; Sag/Cor/Tra), Device Coordinate System (XYZ) and Gradient Coordinate System
    (GCS or PRS; Phase, Readout, Slice) are implemented (so far only rotation, i.e. won't work for offcenter measurementes).

    Examples
    ----------
    ```
    import twixtools
    twix = twixtools.read_twix('meas.dat', parse_geometry=True, parse_data=False)
    x = [1,1,1]
    y = twix[-1]['geometry'].rps_to_xyz() @ x
    ```

    Based on work from Christian Mirkes and Ali Aghaeifar.
    """

    def __init__(self, twix):
        self.from_twix(twix)

    def __str__(self):
        return ("Geometry:\n"
                f"  inplane_rot: {self.inplane_rot}\n"
                f"  normal: {self.normal}\n"
                f"  offset: {self.offset}\n"
                f"  patient_position: {self.patient_position}\n"
                f"  rotmatrix: {self.rotmatrix}\n"
                f"  voxelsize: {self.voxelsize}")

    def from_twix(self, twix):
        if twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 2:
            self.dims = 2
        elif twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 4:
            self.dims = 3
        else:
            self.dims = None

        if len(twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"]) > 1:
            print("WARNING more than one slice. Taking first one..")

        self.fov = [
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dReadoutFOV"]
            * internal_os,
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dPhaseFOV"],
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]["dThickness"],
        ]

        self.resolution = [
            twix["hdr"]["MeasYaps"]["sKSpace"]["lBaseResolution"] * internal_os,
            twix["hdr"]["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
            twix["hdr"]["MeasYaps"]["sKSpace"]["lPartitions"] if self.dims == 3 else 1,
        ]

        self.voxelsize = list(np.array(self.fov) / np.array(self.resolution))

        self.normal = [0, 0, 0]
        if "sNormal" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]:
            for i, d in enumerate(pcs_directions):
                self.normal[i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0][
                    "sNormal"
                ].get(d, self.normal[i])

        self.offset = [0, 0, 0]
        if "sPosition" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0]:
            for i, d in enumerate(pcs_directions):
                self.offset[i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0][
                    "sPosition"
                ].get(d, self.offset[i])

        self.inplane_rot = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][0].get(
            "dInPlaneRot", 0
        )

        if "tPatientPosition" in twix["hdr"]["Meas"]:
            self.patient_position = twix["hdr"]["Meas"].get("tPatientPosition")
        elif "sPatPosition" in twix["hdr"]["Meas"]:
            self.patient_position = twix["hdr"]["Meas"].get("sPatPosition")
        else:
            self.patient_position = None

        self.rotmatrix = self.rps_to_xyz().tolist()

    def get_plane_orientation(self):
        # sanity check if normal vector is unit vector
        norm = np.linalg.norm(self.normal)
        if not abs(1 - norm) < 0.001:
            raise RuntimeError(f"Normal vector is not normal: |x| = {norm}")

        # find main direction of normal vector for first part of rot matrix
        maindir = np.argmax(np.abs(self.normal))
        if maindir == 0:
            init_mat = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]  # @ mat // inplane mat
        elif maindir == 1:
            init_mat = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
        else:
            init_mat = np.eye(3)

        # initialize normal vector direction to which to compute the second part of rotation matrix
        init_normal = np.zeros(3)
        init_normal[maindir] = 1

        # calculate cross product and sine, cosine
        v = np.cross(init_normal, self.normal)
        s = np.linalg.norm(v)
        c = np.dot(init_normal, self.normal)

        if s <= 0.00001:
            # we have cosine 1 or -1, two vectors are (anti-) parallel
            mat = np.matmul(np.eye(3) * c, init_mat)
        else:
            # calculate cross product matrix
            v_x = np.cross(np.eye(3), v)
            # calculate rotation matrix, division should be possible from excluding c = -1 above
            mat = np.eye(3) + v_x + np.divide(np.matmul(v_x, v_x), 1 + c)
            # calculate full rotation matrix
            mat = np.matmul(mat, init_mat)

        return mat

    def get_inplane_rotation(self):
        mat = [
            [-np.sin(self.inplane_rot), np.cos(self.inplane_rot), 0],
            [-np.cos(self.inplane_rot), -np.sin(self.inplane_rot), 0],
            [0, 0, 1],
        ]
        return np.array(mat)

    def prs_to_pcs(self):
        mat = self.get_inplane_rotation()
        mat = self.get_plane_orientation() @ mat
        return mat

    def pcs_to_xyz(self):
        if self.patient_position in pcs_transformations:
            return np.array(pcs_transformations[self.patient_position])
        else:
            raise RuntimeError(f"Unknown patient position: {self.patient_position}")

    def prs_to_xyz(self):
        return self.pcs_to_xyz() @ self.prs_to_pcs()

    def rps_to_xyz(self):
        return self.prs_to_xyz() @ self.rps_to_prs()

    def rps_to_prs(self):
        return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
