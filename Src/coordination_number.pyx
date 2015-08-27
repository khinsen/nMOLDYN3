include 'python.pxi'
include 'numeric.pxi'

cdef extern from "math.h":

    double floor(double x)
    double ceil(double x)
    double sqrt(double x)

cdef inline double round(double r):
    return floor(r + 0.5) if (r > 0.0) else ceil(r - 0.5)

def coordination_number(array_type configuration, array_type directCell, array_type reverseCell,\
                        array_type gIndexes, int mol, array_type indexes, array_type molecules, array_type elements,\
                        array_type cnIntra, array_type cnInter, array_type scalec, array_type groupcenter,\
                        double rmin, double dr):

    cdef int *atomindex, *molindex, *elementindex, *groupindex
    cdef double *config, *cell, *rcell, *scaleconfig, *center
    cdef float* coordNumberIntra, *coordNumberInter

    cdef double x, y, z, sdx, sdy, sdz, rx, ry, rz, r

    cdef int i, ind, bin, index, dim, nbins

    # Checks the dimensions of the input arrays.
    assert configuration.nd == 2
    assert directCell.nd == 1
    assert reverseCell.nd == 1
    assert gIndexes.nd == 1
    assert indexes.nd == 1
    assert molecules.nd == 1    
    assert elements.nd == 1
    assert cnIntra.nd == 2
    assert cnInter.nd == 2
    assert scalec.nd == 1
    assert groupcenter.nd == 1
    
    # Checks the types of the input arrays.
    assert configuration.descr.type_num == PyArray_DOUBLE
    assert directCell.descr.type_num    == PyArray_DOUBLE
    assert reverseCell.descr.type_num   == PyArray_DOUBLE
    assert scalec.descr.type_num        == PyArray_DOUBLE
    assert groupcenter.descr.type_num   == PyArray_DOUBLE
    
    assert (gIndexes.descr.type_num == PyArray_INT or gIndexes.descr.type_num == PyArray_LONG) \
    and gIndexes.descr.elsize == sizeof(int)

    assert (indexes.descr.type_num == PyArray_INT or indexes.descr.type_num == PyArray_LONG) \
    and indexes.descr.elsize == sizeof(int)
    
    assert (molecules.descr.type_num == PyArray_INT or molecules.descr.type_num == PyArray_LONG) \
    and molecules.descr.elsize == sizeof(int)
    
    assert (elements.descr.type_num == PyArray_INT or elements.descr.type_num == PyArray_LONG) \
    and elements.descr.elsize == sizeof(int)

    assert cnIntra.descr.type_num == PyArray_FLOAT
    assert cnInter.descr.type_num == PyArray_FLOAT

    assert cnInter.dimensions[0] == cnIntra.dimensions[0]
    assert cnInter.dimensions[1] == cnIntra.dimensions[1]

    nbins = cnInter.dimensions[1]

    config = <double *>configuration.data
    cell = <double *>directCell.data
    rcell = <double *>reverseCell.data
    groupindex = <int *>gIndexes.data
    atomindex = <int *>indexes.data
    molindex = <int *>molecules.data
    elementindex = <int *>elements.data
    coordNumberIntra = <float *>cnIntra.data
    coordNumberInter = <float *>cnInter.data
    scaleconfig = <double *> scalec.data
    center = <double *> groupcenter.data
    
    center[0] = 0.0
    center[1] = 0.0
    center[2] = 0.0

    for 0 <= i < gIndexes.dimensions[0]:

        ind = 3*groupindex[i]

        x = config[ind]
        y = config[ind+1]
        z = config[ind+2]

        center[0] += (x*rcell[0] + y*rcell[3] + z*rcell[6])
        center[1] += (x*rcell[1] + y*rcell[4] + z*rcell[7])
        center[2] += (x*rcell[2] + y*rcell[5] + z*rcell[8])

    center[0] /= <double>gIndexes.dimensions[0]
    center[1] /= <double>gIndexes.dimensions[0]
    center[2] /= <double>gIndexes.dimensions[0]    
                   
    for 0 <= i < indexes.dimensions[0]:

        ind = 3*atomindex[i]

        x = config[ind]
        y = config[ind+1]
        z = config[ind+2]

        ind = 3*i

        scaleconfig[ind]   = x*rcell[0] + y*rcell[3] + z*rcell[6]
        scaleconfig[ind+1] = x*rcell[1] + y*rcell[4] + z*rcell[7]
        scaleconfig[ind+2] = x*rcell[2] + y*rcell[5] + z*rcell[8]

    dim  = cnInter.dimensions[1]
    
    for 0 <= i < indexes.dimensions[0]:
        
        ind = 3*i

        sdx = scaleconfig[ind]   - center[0]
        sdy = scaleconfig[ind+1] - center[1]
        sdz = scaleconfig[ind+2] - center[2]

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
                    
        index = elementindex[i]*dim + bin
                
        if mol == molindex[i]:
            coordNumberIntra[index] += 1.0
        else:
            coordNumberInter[index] += 1.0
            
