include 'python.pxi'
include 'numeric.pxi'

cdef extern from "math.h":

    double floor(double x)
    double ceil(double x)
    double sqrt(double x)

cdef inline double round(double r):
    return floor(r + 0.5) if (r > 0.0) else ceil(r - 0.5)

def distance_histogram(array_type configuration, array_type directCell, array_type reverseCell,\
                       array_type indexes, array_type molecules, array_type elements,\
                       array_type histintra, array_type histinter,\
                       array_type scalec, double rmin, double dr):

    # This computes the intra and intermolecular distances histogram.
    # The algorithm is a Pyrex adaptation of the FORTRAN implementation 
    # made by Miguel Angel Gonzalez (Institut Laue Langevin).

    cdef int *atomindex, *molindex, *elementindex
    cdef double *config, *cell, *rcell, *scaleconfig
    cdef float *hintra, *hinter

    cdef double x, y, z, sdx, sdy, sdz, rx, ry, rz, r

    cdef int i, j, ind, indi, indj, bin, hindex, dimjk, dimk, nbins

    # Checks the dimensions of the input arrays.
    assert configuration.nd == 2
    assert directCell.nd == 1
    assert reverseCell.nd == 1
    assert indexes.nd == 1
    assert molecules.nd == 1
    assert elements.nd == 1
    assert histintra.nd == 3
    assert histinter.nd == 3
    assert scalec.nd == 1

    # Checks the types of the input arrays.
    assert configuration.descr.type_num == PyArray_DOUBLE
    assert directCell.descr.type_num    == PyArray_DOUBLE
    assert reverseCell.descr.type_num   == PyArray_DOUBLE
    assert scalec.descr.type_num        == PyArray_DOUBLE

    assert (indexes.descr.type_num == PyArray_INT or indexes.descr.type_num == PyArray_LONG) \
    and indexes.descr.elsize == sizeof(int)
    
    assert (molecules.descr.type_num == PyArray_INT or molecules.descr.type_num == PyArray_LONG) \
    and molecules.descr.elsize == sizeof(int)
    
    assert (elements.descr.type_num == PyArray_INT or elements.descr.type_num == PyArray_LONG) \
    and elements.descr.elsize == sizeof(int)
    
    assert histintra.descr.type_num == PyArray_FLOAT
    assert histinter.descr.type_num == PyArray_FLOAT

    assert histinter.dimensions[0] == histintra.dimensions[0]
    assert histinter.dimensions[1] == histintra.dimensions[1]
    assert histinter.dimensions[2] == histintra.dimensions[2]

    nbins = histinter.dimensions[2]

    config = <double *>configuration.data
    cell = <double *>directCell.data
    rcell = <double *>reverseCell.data
    atomindex = <int *>indexes.data
    molindex = <int *>molecules.data
    elementindex = <int *>elements.data
    hintra = <float *>histintra.data
    hinter = <float *>histinter.data
    scaleconfig = <double *> scalec.data
    
    for 0 <= i < indexes.dimensions[0]:

        ind = 3*atomindex[i]

        x = config[ind]
        y = config[ind+1]
        z = config[ind+2]

        ind = 3*i

        scaleconfig[ind]   = x*rcell[0] + y*rcell[3] + z*rcell[6]
        scaleconfig[ind+1] = x*rcell[1] + y*rcell[4] + z*rcell[7]
        scaleconfig[ind+2] = x*rcell[2] + y*rcell[5] + z*rcell[8]
            
    dimjk = histinter.dimensions[1]*histinter.dimensions[2] 
    dimk  = histinter.dimensions[2]

    for 0 <= i < indexes.dimensions[0] - 1:

        indi = 3*i

        sx = scaleconfig[indi]
        sy = scaleconfig[indi+1]
        sz = scaleconfig[indi+2]

        for i + 1 <= j < indexes.dimensions[0]:

            indj = 3*j

            sdx = scaleconfig[indj]   - sx
            sdy = scaleconfig[indj+1] - sy
            sdz = scaleconfig[indj+2] - sz

            sdx -= round(sdx)
            sdy -= round(sdy)
            sdz -= round(sdz)
            
            rx = sdx*cell[0] + sdy*cell[3] + sdz*cell[6]
            ry = sdx*cell[1] + sdy*cell[4] + sdz*cell[7]
            rz = sdx*cell[2] + sdy*cell[5] + sdz*cell[8]

            r = sqrt(rx*rx + ry*ry + rz*rz)

            bin = <int>((r-rmin)/dr)

            if bin < 0 or bin >= nbins:
                continue

            hindex = elementindex[i]*dimjk + elementindex[j]*dimk + bin

            if molindex[i] == molindex[j]:
                hintra[hindex] += 1.0
            else:
                hinter[hindex] += 1.0
