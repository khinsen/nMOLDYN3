include 'python.pxi'
include 'numeric.pxi'

cdef extern from "math.h":

    double round(double x)
    double sqrt(double x)
    double exp(double x)
    double tanh(double x)

def order_parameter_pyrex(array_type configuration,\
                          int hIndex,\
                          int oIndex,\
                          array_type directCell,\
                          array_type reverseCell,\
                          array_type indexes,\
                          array_type scalec,\
                          array_type scaleh,\
                          array_type scaleo,\
                          double conversionFactor):

    cdef int *atomindex
    cdef double *config, *cell, *rcell, *scaleconfig, *scalehposition, *scaleoposition

    cdef double x, y, z, sdx, sdy, sdz, rx, ry, rz, r, s2

    cdef int i, ind

    # Checks the dimension of the input arrays.
    assert configuration.nd == 2
    assert directCell.nd == 1
    assert reverseCell.nd == 1
    assert indexes.nd == 1
    assert scalec.nd == 1
    assert scaleh.nd == 1
    assert scaleo.nd == 1

    # Checks the types of the input arrays.
    assert configuration.descr.type_num  == PyArray_DOUBLE
    assert directCell.descr.type_num     == PyArray_DOUBLE
    assert reverseCell.descr.type_num    == PyArray_DOUBLE
    assert scalec.descr.type_num         == PyArray_DOUBLE
    assert scaleh.descr.type_num         == PyArray_DOUBLE
    assert scaleo.descr.type_num         == PyArray_DOUBLE

    assert (indexes.descr.type_num == PyArray_INT or indexes.descr.type_num == PyArray_LONG) \
    and indexes.descr.elsize == sizeof(int)
        
    config = <double *>configuration.data
    cell = <double *>directCell.data
    rcell = <double *>reverseCell.data
    atomindex = <int *>indexes.data
    scaleconfig = <double *> scalec.data
    scalehposition = <double *> scaleh.data
    scaleoposition = <double *> scaleo.data
        
    scalehposition[0] = config[3*hIndex]*rcell[0] + config[3*hIndex+1]*rcell[3] + config[3*hIndex+2]*rcell[6]
    scalehposition[1] = config[3*hIndex]*rcell[1] + config[3*hIndex+1]*rcell[4] + config[3*hIndex+2]*rcell[7]
    scalehposition[2] = config[3*hIndex]*rcell[2] + config[3*hIndex+1]*rcell[5] + config[3*hIndex+2]*rcell[8]

    scaleoposition[0] = config[3*oIndex]*rcell[0] + config[3*oIndex+1]*rcell[3] + config[3*oIndex+2]*rcell[6]
    scaleoposition[1] = config[3*oIndex]*rcell[1] + config[3*oIndex+1]*rcell[4] + config[3*oIndex+2]*rcell[7]
    scaleoposition[2] = config[3*oIndex]*rcell[2] + config[3*oIndex+1]*rcell[5] + config[3*oIndex+2]*rcell[8]
        
    for 0 <= i < indexes.dimensions[0]:

        ind = 3*atomindex[i]

        x = config[ind]
        y = config[ind+1]
        z = config[ind+2]

        ind = 3*i

        scaleconfig[ind]   = x*rcell[0] + y*rcell[3] + z*rcell[6]
        scaleconfig[ind+1] = x*rcell[1] + y*rcell[4] + z*rcell[7]
        scaleconfig[ind+2] = x*rcell[2] + y*rcell[5] + z*rcell[8]
            
    s2 = 0.0
    for 0 <= i < indexes.dimensions[0]:

        ind = 3*i

        sdx = scaleconfig[ind]   - scalehposition[0]
        sdy = scaleconfig[ind+1] - scalehposition[1]
        sdz = scaleconfig[ind+2] - scalehposition[2]

        sdx -= round(sdx)
        sdy -= round(sdy)
        sdz -= round(sdz)
            
        rx = sdx*cell[0] + sdy*cell[3] + sdz*cell[6]
        ry = sdx*cell[1] + sdy*cell[4] + sdz*cell[7]
        rz = sdx*cell[2] + sdy*cell[5] + sdz*cell[8]

        r = sqrt(rx*rx + ry*ry + rz*rz)
        
        r /= conversionFactor

        s2 += 0.8*exp(-r)
            
        sdx = scaleconfig[ind]   - scaleoposition[0]
        sdy = scaleconfig[ind+1] - scaleoposition[1]
        sdz = scaleconfig[ind+2] - scaleoposition[2]

        sdx -= round(sdx)
        sdy -= round(sdy)
        sdz -= round(sdz)
            
        rx = sdx*cell[0] + sdy*cell[3] + sdz*cell[6]
        ry = sdx*cell[1] + sdy*cell[4] + sdz*cell[7]
        rz = sdx*cell[2] + sdy*cell[5] + sdz*cell[8]

        r = sqrt(rx*rx + ry*ry + rz*rz)

        r /= conversionFactor

        s2 += 1.0*exp(-r)
        
    s2 = tanh(2.656*s2) - 0.1
        
    return s2

