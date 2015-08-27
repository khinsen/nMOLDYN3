include 'python.pxi'
include 'numeric.pxi'

cdef extern from "math.h":

    double round(double x)
    double sin(double x)
    double sqrt(double x)

def smoothed_static_coherent_structure_factor_pyrex(array_type configuration,\
                                                    array_type directCell,\
                                                    array_type reverseCell,\
                                                    array_type qvalues,\
                                                    array_type indexes,\
                                                    array_type elements,\
                                                    array_type ssftemp,\
                                                    array_type scalec):

    cdef int *atomindex, *elementindex
    cdef double *config, *cell, *rcell, *scaleconfig, *qval, *ssf

    cdef double x, y, z, sdx, sdy, sdz, rx, ry, rz, r, qr

    cdef int i, j, ind, indi, indj, q, dimjk, dimk, index
    
    # Checks the dimension of the input arrays.
    assert configuration.nd == 2
    assert directCell.nd == 1
    assert reverseCell.nd == 1
    assert qvalues.nd == 1
    assert indexes.nd == 1
    assert elements.nd == 1
    assert ssftemp.nd == 3
    assert scalec.nd == 1

    # Checks the types of the input arrays.
    assert configuration.descr.type_num == PyArray_DOUBLE
    assert directCell.descr.type_num    == PyArray_DOUBLE
    assert reverseCell.descr.type_num   == PyArray_DOUBLE
    assert qvalues.descr.type_num       == PyArray_DOUBLE
    assert scalec.descr.type_num        == PyArray_DOUBLE

    assert (indexes.descr.type_num == PyArray_INT or indexes.descr.type_num == PyArray_LONG) \
    and indexes.descr.elsize == sizeof(int)
        
    assert (elements.descr.type_num == PyArray_INT or elements.descr.type_num == PyArray_LONG) \
    and elements.descr.elsize == sizeof(int)
    
    config = <double *>configuration.data
    cell = <double *>directCell.data
    rcell = <double *>reverseCell.data
    qval = <double *>qvalues.data
    atomindex = <int *>indexes.data
    elementindex = <int *>elements.data
    ssf = <double *>ssftemp.data
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

    dimjk = ssftemp.dimensions[1]*ssftemp.dimensions[2] 
    dimk  = ssftemp.dimensions[2]

    for 0 <= i < indexes.dimensions[0] - 1:

        indi = 3*i

        x = scaleconfig[indi]
        y = scaleconfig[indi+1]
        z = scaleconfig[indi+2]

        for i + 1 <= j < indexes.dimensions[0]:
                
            indj = 3*j

            sdx = scaleconfig[indj]   - x
            sdy = scaleconfig[indj+1] - y
            sdz = scaleconfig[indj+2] - z

            sdx -= round(sdx)
            sdy -= round(sdy)
            sdz -= round(sdz)

            rx = sdx*cell[0] + sdy*cell[3] + sdz*cell[6]
            ry = sdx*cell[1] + sdy*cell[4] + sdz*cell[7]
            rz = sdx*cell[2] + sdy*cell[5] + sdz*cell[8]

            r = sqrt(rx*rx + ry*ry + rz*rz)
            
            index = elementindex[i]*dimjk + elementindex[j]*dimk
            
            for 0 <= k < qvalues.dimensions[0]:
                        
                qr = qval[k]*r
                
                if qr == 0:
                    ssf[index] += 1.0
                
                else:
                    ssf[index] += sin(qr)/qr

                index += 1
