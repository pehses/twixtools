#!/usr/bin/python
import numpy as np
import warnings
internal_os = 2
pcs_directions = ["dSag", "dCor", "dTra"]

# p. 418 - pcs to dcs
pcs_transformations = {
    "HFS": [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    "HFP": [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    "FFS": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
}


# scalar first convention - see MDB.
def quat_to_rotmat(scalar, i, j, k):
    quat = np.array([scalar, i, j, k])

    norm = np.linalg.norm(quat)
    if abs(1 - norm) < 1e-6:
        warnings.warn(f"Quaternion is not normalized (norm = {norm})", UserWarning)

    r = scalar

    mat = np.array([
            [ 1 - 2 * (j**2 + k**2), 2 * (i * j - k * r)  , 2 * (i * k + j * r)   ],
            [ 2 * (i * j + k * r)  , 1 - 2 * (i**2 + k**2), 2 * (j * k - i * r)   ],
            [ 2 * (i * k - j * r)  , 2 * (j * k + i * r)  , 1 - 2 * (i**2 + j**2) ]])
    return mat


def parse_slice_order(twix):
    order = None
    if '-' != twix['hdr']['Config']['chronSliceIndices'][0]:
        order = []
        for x in twix['hdr']['Config']['chronSliceIndices']:
            if len(order) == int(twix['hdr']['MeasYaps']['sSliceArray']['lSize']):
                break
            if x == ' ':
                continue
            val = int(x)
            order.append(val)
    return order


def prs2sct_mdb(twix, sliceno):
    """Extract orientation matrix from mdb"""

    # match chronological and normal slice order:
    order = parse_slice_order(twix)
    if order:
        lookup = { x: i for i, x in enumerate(order) }
        original_index = sliceno
        sliceno = lookup[sliceno]

    # find first mdb which belongs to the slice:
    index = -1
    for i, m in enumerate(twix['mdb']):
        if m.mdh.Counter.Sli == sliceno:
            index = i
            break

    if -1 == index:
        print(f"No MDB found for slice with chron. index {original_index}, index {sliceno}.")
        raise RuntimeError("Geom-MDB-Not-Found")

    # calc rot matrix:
    mat = quat_to_rotmat(*twix['mdb'][index].mdh.SliceData.Quaternion)

    # readout and pe are flipped:
    mat[:, :2] *= -1

    return mat


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
    """

    @staticmethod
    def create_for_all_slices(twix):
        slices = len(twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"])
        return [Geometry(twix, n_slice=i) for i in range(slices)]

    def __init__(self, twix, n_slice=None):
        self.from_twix(twix, n_slice)

    def __str__(self):
        return ("Geometry:\n"
                f"  inplane_rot: {self.inplane_rot}\n"
                f"  normal: {self.normal}\n"
                f"  offset: {self.offset}\n"
                f"  patient_position: {self.patient_position}\n"
                f"  voxelsize: {self.voxelsize}")

    def from_twix(self, twix, n_slice=None):
        if twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 2:
            self.dims = 2
        elif twix["hdr"]["MeasYaps"]["sKSpace"]["ucDimension"] == 4:
            self.dims = 3
        else:
            self.dims = None

        if n_slice is None and len(twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"]) > 1:
            print("WARNING: Geometry calculations are valid only for the first slice in this multi-slice acquisition.")
            n_slice = 0

        self.fov = [
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice]["dReadoutFOV"]
            * internal_os,
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice]["dPhaseFOV"],
            twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice]["dThickness"],
        ]

        self.resolution = [
            twix["hdr"]["MeasYaps"]["sKSpace"]["lBaseResolution"] * internal_os,
            twix["hdr"]["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
            twix["hdr"]["MeasYaps"]["sKSpace"]["lPartitions"] if self.dims == 3 else 1,
        ]

        self.voxelsize = list(np.array(self.fov) / np.array(self.resolution))

        self.normal = [0, 0, 0]
        if "sNormal" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice]:
            for i, d in enumerate(pcs_directions):
                self.normal[i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice][
                    "sNormal"
                ].get(d, self.normal[i])

        if not 'mdb' in twix or 0 == len(twix['mdb']):
            print("Trying to create slice geometry, but no data was found \
                    (i.e. parse_geometry AND NOT parse_data.) - this is not possible.\n\
                    Please either set parse_geometry = False or parse_data = True")

            raise RuntimeError("MDBs-Needed-For-Slice-Geometry")

        self._prs_to_pcs_mat = prs2sct_mdb(twix, n_slice)

        self.offset = [0, 0, 0]
        if "sPosition" in twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice]:
            for i, d in enumerate(pcs_directions):
                self.offset[i] = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice][
                    "sPosition"
                ].get(d, self.offset[i])

        self.inplane_rot = twix["hdr"]["MeasYaps"]["sSliceArray"]["asSlice"][n_slice].get(
            "dInPlaneRot", 0
        )

        if "tPatientPosition" in twix["hdr"]["Meas"]:
            self.patient_position = twix["hdr"]["Meas"].get("tPatientPosition")
        elif "sPatPosition" in twix["hdr"]["Meas"]:
            self.patient_position = twix["hdr"]["Meas"].get("sPatPosition")
        else:
            self.patient_position = None

    def prs_to_pcs(self):
        return self._prs_to_pcs_mat

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
