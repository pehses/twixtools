#!/usr/bin/python
import numpy as np

def apply_transform(g, data):
    import scipy.ndimage
    data = resample(g, data)
    a = _augment_matrix_dims(np.array(g['matrix']), len(data.shape))
    A = np.array(np.linalg.inv(a))
    output_shape,offset = get_newshape(data.shape, a)
    return scipy.ndimage.affine_transform(data, A, offset=offset, output_shape=output_shape)

def get_newshape(shape, mat):
    edges=get_edges(shape)

    for i,edge in enumerate(edges):
        edges[i] = mat @ edge
    #get bottom left edge.
    bottom = np.amin(edges, axis=0)
    shifted_edges = np.array([edge - bottom for edge in edges])
    # new top right edge
    output_shape = np.amax(shifted_edges, axis=0)
    offset = np.array(shape)/2 - np.linalg.inv(mat) @ (output_shape/2)

    output_shape = tuple([int(x) for x in output_shape])
    print(output_shape, offset)

    return output_shape, offset

def get_edges(shape):
    edges = []
    # iterate over all possible selections of dimensions
    for bitmask in range(2**len(shape)):
        edge = np.zeros(len(shape))
        for dim in range(len(shape)):
            if bitmask & 1<<dim:
                edge[dim]=shape[dim]
        edges.append(edge)
    return edges

def resample(g, data):
    import scipy.interpolate
    def get_isotropic_dims(g):
        ref_vs = np.mean(g['voxelsize'])
        dims = np.copy(g['resolution'])
        for i,vs in enumerate(g['voxelsize']):
            if abs(vs - ref_vs) > 5 / dims[i] * vs:
                dims[i] = dims[i] * vs / ref_vs
        return dims

    newdims = get_isotropic_dims(g)
    maxdist = max(np.abs(g['resolution']) - np.array(newdims))
    if maxdist > 5:
        print(f"""
        Old dims: {g['resolution']}
        Old voxelsize: {g['voxelsize']};
        New dims: {newdims};
        New voxelsize: {g['fov']/newdims}
        """)
        g['iso_dims'] = newdims.tolist()
        g['iso_voxelsize'] = (g['fov']/newdims).tolist()

        dims = data.shape
        coor = [ (np.arange(s) - (s-1)/2)/s for s in dims]
        coor_neu = [ (np.arange(s) - (s-1)/2)/s for s in newdims]
        target_grid = np.stack(np.meshgrid(*coor_neu, indexing='ij'), axis=-1)
        return scipy.interpolate.interpn(tuple(coor), data, target_grid, bounds_error=False, fill_value=None, method='linear')
    else:
        return data

def _augment_matrix_dims(a, N):
    e = np.eye(N)
    e[:a.shape[0],:a.shape[1]] = a[:,:]
    return e


def create_transform_file(transform, raw):
    import twixtools
    from twixtools import geometry

    twix = twixtools.read_twix(raw, parse_data=False)
    g = geometry.get_geometry(twix[-1])

    with open(transform, 'w') as f:
        json.dump(g, f, indent=4)

def apply_transform_file(transform, infile, outfile):
    from twixtools.contrib import cfl
    data = cfl.readcfl(infile)
    with open(transform,'r') as f:
        g = json.load(f)

    data = apply_transform(g, data)
    cfl.writecfl(outfile, data)

if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("transform", help="JSON File for transformation")
    parser.add_argument("infile", help="Either .dat file or cfl input file to be transformed")
    parser.add_argument("outfile", nargs="?", default=None, help="Transformed output (optional)")
    args = parser.parse_args()

    if args.outfile:
        apply_transform_file(args.transform, args.infile, args.outfile)
    else:
        create_transform_file(args.transform, args.infile)
