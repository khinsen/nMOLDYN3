include 'python.pxi'
include 'numeric.pxi'

ctypedef int size_t

cdef extern from "stdlib.h":
    void *malloc(size_t size)
    int free(void*)
    int sizeof()

cdef extern from "math.h":

    double floor(double x)
    double ceil(double x)
    double sqrt(double x)
    double cos(double x)
    double acos(double x)
    double fmin(double x, double y)
    double fmax(double x, double y)

cdef inline double round(double r):
    return floor(r + 0.5) if (r > 0.0) else ceil(r - 0.5)
            
# H-Bond detection algorithm in case of a periodic universe.
def hbond_detection(array_type configuration, array_type directCell, array_type reverseCell,\
                    array_type acc_indexes, array_type don_indexes, array_type h_indexes,\
                    array_type hbond_matrix, double dis_min, double dis_max, double ang_min, double ang_max):

    cdef int *acc_idx, *don_idx, *h_idx
    
    cdef double *config, *cell, *rcell, *scfg_acc, *scfg_don, *scfg_h
    
    cdef float *hbond_mat
    
    cdef int i, j, ind
    
    cdef double rx, ry, rz, r1, cosa, angle, \
                x_acc, y_acc, z_acc, s_x_acc, s_y_acc, s_z_acc, \
                x_don, y_don, z_don, s_x_don, s_y_don, s_z_don, \
                x_h, y_h, z_h, s_x_h, s_y_h, s_z_h, \
                dx, dy, dz, dx1, dy1, dz1
    
    # Checks the dimensions of the input arrays.
    assert configuration.nd == 2
    assert directCell.nd    == 1
    assert reverseCell.nd   == 1
    assert acc_indexes.nd   == 1
    assert don_indexes.nd   == 1
    assert h_indexes.nd     == 1
    assert don_indexes.dimensions[0] == h_indexes.dimensions[0]
    assert hbond_matrix.nd  == 2

    # Checks the types of the input arrays.
    assert configuration.descr.type_num == PyArray_DOUBLE
    assert directCell.descr.type_num    == PyArray_DOUBLE
    assert reverseCell.descr.type_num   == PyArray_DOUBLE
    assert (acc_indexes.descr.type_num == PyArray_INT or acc_indexes.descr.type_num == PyArray_LONG) \
    and acc_indexes.descr.elsize == sizeof(int)
    assert (don_indexes.descr.type_num == PyArray_INT or don_indexes.descr.type_num == PyArray_LONG) \
    and don_indexes.descr.elsize == sizeof(int)
    assert (h_indexes.descr.type_num == PyArray_INT or h_indexes.descr.type_num == PyArray_LONG) \
    and h_indexes.descr.elsize == sizeof(int)        
    assert hbond_matrix.descr.type_num == PyArray_FLOAT

    # The input array parameters are assigned to local C array variables.
    config     = <double *>configuration.data
    cell       = <double *>directCell.data
    rcell      = <double *>reverseCell.data
    acc_idx    = <int *>acc_indexes.data
    don_idx    = <int *>don_indexes.data
    h_idx      = <int *>h_indexes.data
    hbond_mat  = <float *>hbond_matrix.data
    
    # The C array that will store the scaled coordinates of acceptors.
    scfg_acc = <double *> malloc(3*acc_indexes.dimensions[0]*sizeof(double))
    
    # The C array that will store the scaled coordinates of donors.
    scfg_don = <double *> malloc(3*don_indexes.dimensions[0]*sizeof(double))
    
    # The C array that will store the scaled coordinates of hydrogens.
    scfg_h   = <double *> malloc(3*h_indexes.dimensions[0]*sizeof(double))
        
    # Loop over the acceptor indexes to define the corresponding scaled configuration.
    for 0 <= i < acc_indexes.dimensions[0]:

        ind = 3*acc_idx[i]

        # The acceptor real coordinates.
        x_acc = config[ind]
        y_acc = config[ind+1]
        z_acc = config[ind+2]

        ind = 3*i

        # The acceptor scaled coordinates are computed.
        scfg_acc[ind]   = x_acc*rcell[0] + y_acc*rcell[3] + z_acc*rcell[6]
        scfg_acc[ind+1] = x_acc*rcell[1] + y_acc*rcell[4] + z_acc*rcell[7]
        scfg_acc[ind+2] = x_acc*rcell[2] + y_acc*rcell[5] + z_acc*rcell[8]

    # Loop over the donor and hydrogen indexes to define the corresponding scaled configuration.
    for 0 <= i < don_indexes.dimensions[0]:

        ind = 3*don_idx[i]

        # The donor real coordinates.
        x_don = config[ind]
        y_don = config[ind+1]
        z_don = config[ind+2]

        ind = 3*i

        # The donor scaled coordinates are computed.
        scfg_don[ind]   = x_don*rcell[0] + y_don*rcell[3] + z_don*rcell[6]
        scfg_don[ind+1] = x_don*rcell[1] + y_don*rcell[4] + z_don*rcell[7]
        scfg_don[ind+2] = x_don*rcell[2] + y_don*rcell[5] + z_don*rcell[8]
        
        ind = 3*h_idx[i]

        # The hydrogen real coordinates.
        x_h = config[ind]
        y_h = config[ind+1]
        z_h = config[ind+2]

        ind = 3*i

        # The hydrogen scaled coordinates are computed.
        scfg_h[ind]   = x_h*rcell[0] + y_h*rcell[3] + z_h*rcell[6]
        scfg_h[ind+1] = x_h*rcell[1] + y_h*rcell[4] + z_h*rcell[7]
        scfg_h[ind+2] = x_h*rcell[2] + y_h*rcell[5] + z_h*rcell[8]
        
    # Loop over the acceptor atoms.
    for 0 <= i < acc_indexes.dimensions[0]:
    
        ind = 3*acc_idx[i]

        # The acceptor real coordinates.
        x_acc = config[ind]
        y_acc = config[ind+1]
        z_acc = config[ind+2]

        ind = 3*i

        # The acceptor scaled coordinates.
        s_x_acc = scfg_acc[ind]
        s_y_acc = scfg_acc[ind+1]
        s_z_acc = scfg_acc[ind+2]
        
        for 0 <= j < don_indexes.dimensions[0]:
        
            ind = 3*j
            
            # The donor scaled real coordinates.
            s_x_don = scfg_don[ind]
            s_y_don = scfg_don[ind+1]
            s_z_don = scfg_don[ind+2]
            
            # The components of the Acc-Don vector in scaled coordinates.
            dx = s_x_don - s_x_acc
            dy = s_y_don - s_y_acc
            dz = s_z_don - s_z_acc
            
            # The PBC are applied.
            dx -= round(dx)
            dy -= round(dy)
            dz -= round(dz)
            
            # The Acc-Don minimum distance in real coordinates is computed.
            rx = dx*cell[0] + dy*cell[3] + dz*cell[6]
            ry = dx*cell[1] + dy*cell[4] + dz*cell[7]
            rz = dx*cell[2] + dy*cell[5] + dz*cell[8]
            r = sqrt(rx*rx + ry*ry + rz*rz)
            
            # If it is not in the H-Bond accepted interval, do not go futher.
            if r < dis_min or r > dis_max: continue
                        
            # The new scaled coordinates for the current donor atom are computed.
            s_x_don = s_x_acc + dx
            s_y_don = s_y_acc + dy
            s_z_don = s_z_acc + dz
            
            # The donor 'unfolded' real coordinates are computed.
            x_don = s_x_don*cell[0] + s_y_don*cell[3] + s_z_don*cell[6]
            y_don = s_x_don*cell[1] + s_y_don*cell[4] + s_z_don*cell[7]
            z_don = s_x_don*cell[2] + s_y_don*cell[5] + s_z_don*cell[8]
            
            # The hydrogen scaled coordinates.                                                
            s_x_h = scfg_h[ind]
            s_y_h = scfg_h[ind+1]
            s_z_h = scfg_h[ind+2]
            
            # The components of the Acc-H vector in scaled coordinates.
            dx = s_x_h - s_x_acc
            dy = s_y_h - s_y_acc
            dz = s_z_h - s_z_acc
            
            # The PBC are applied.
            s_x_h -= round(dx)
            s_y_h -= round(dy)
            s_z_h -= round(dz)
                                    
            # The hydrogen 'unfolded' real coordinates are computed.
            x_h = s_x_h*cell[0] + s_y_h*cell[3] + s_z_h*cell[6]
            y_h = s_x_h*cell[1] + s_y_h*cell[4] + s_z_h*cell[7]
            z_h = s_x_h*cell[2] + s_y_h*cell[5] + s_z_h*cell[8]
                
            # The H-Acc vector components.
            dx = x_acc - x_h
            dy = y_acc - y_h
            dz = z_acc - z_h
            
            # The H-Acc vector squared norm.
            r  = dx*dx + dy*dy + dz*dz                        
            
            # The H-Don vector components.
            dx1 = x_don - x_h
            dy1 = y_don - y_h
            dz1 = z_don - z_h
            
            # The H-Don vector squared norm.            
            r1 = dx1*dx1 + dy1*dy1 + dz1*dz1
                       
            # The (H_Acc,H-Don) cosinus is computed.
            cosa = (dx*dx1) + (dy*dy1) + (dz*dz1)
            cosa /=sqrt(r*r1)                        
            cosa = fmax(-1.0,fmin(1.,cosa))
            
            # The angle corresponding to the computed cosinus.
            angle = acos(cosa)
                        
            # If the angle is not in the H-Bond accepted interval, do not go futher.
            if angle < ang_min or angle > ang_max: continue
                                    
            ind = j + i*hbond_matrix.dimensions[1]
                        
            # The H bond matrix is updated with the new H bond found.
            hbond_mat[ind] = 1.0
                                                
    free(scfg_don)
    free(scfg_acc)
    free(scfg_h)
                        
# H-Bond detection algorithm in case of an infinite universe.
def hbond_detection_no_PBC(array_type configuration,\
                           array_type acc_indexes,\
                           array_type don_indexes,\
                           array_type h_indexes,\
                           array_type hbond_matrix,\
                           double dis_min,\
                           double dis_max,\
                           double ang_min,\
                           double ang_max):

    cdef int *acc_idx, *don_idx, *h_idx
    
    cdef double *config
    
    cdef float *hbond_mat
    
    cdef int i, j, ind
    
    cdef double rx, ry, rz, r1, cosa, angle, x_acc, y_acc, z_acc, x_don, y_don, z_don, x_h, y_h, z_h
    
    # Checks the dimensions of the input arrays.
    assert configuration.nd == 2
    assert acc_indexes.nd   == 1
    assert don_indexes.nd   == 1
    assert h_indexes.nd     == 1
    assert hbond_matrix.nd  == 2

    # Checks the types of the input arrays.
    assert configuration.descr.type_num == PyArray_DOUBLE
    assert (acc_indexes.descr.type_num == PyArray_INT or acc_indexes.descr.type_num == PyArray_LONG) \
    and acc_indexes.descr.elsize == sizeof(int)    
    assert (don_indexes.descr.type_num == PyArray_INT or don_indexes.descr.type_num == PyArray_LONG) \
    and don_indexes.descr.elsize == sizeof(int)
    assert (h_indexes.descr.type_num == PyArray_INT or h_indexes.descr.type_num == PyArray_LONG) \
    and h_indexes.descr.elsize == sizeof(int)
    assert hbond_matrix.descr.type_num == PyArray_FLOAT

    # The input array parameters are assigned to local C array variables.
    config     = <double *>configuration.data
    acc_idx    = <int *>acc_indexes.data
    don_idx    = <int *>don_indexes.data
    h_idx      = <int *>h_indexes.data
    hbond_mat  = <float *>hbond_matrix.data
                    
    # Loop over the acceptor atoms.
    for 0 <= i < acc_indexes.dimensions[0]:
    
        ind = 3*acc_idx[i]

        # The acceptor real coordinates.
        x_acc = config[ind]
        y_acc = config[ind+1]
        z_acc = config[ind+2]

        # Loop over the donor atoms.
        for 0 <= j < don_indexes.dimensions[0]:
        
            ind = 3*don_idx[j]
            
            # The donor real coordinates.
            x_don = config[ind]
            y_don = config[ind+1]
            z_don = config[ind+2]
            
            # The components of the Acc-Don vector in scaled coordinates.
            dx = x_don - x_acc
            dy = y_don - y_acc
            dz = z_don - z_acc
            
            # The norm of the Acc-Don vector.
            r = sqrt(dx*dx + dy*dy + dz*dz)
            
            # If it is not in the H-Bond accepted interval, do not go futher.
            if r < dis_min or r > dis_max: continue
                        
            ind = 3*h_idx[j]
                        
            # The hydrogen real coordinates.
            x_h = config[ind]
            y_h = config[ind+1]
            z_h = config[ind+2]
                                                                                       
            # The H-Acc vector components in real coordinates.
            dx = x_acc - x_h
            dy = y_acc - y_h
            dz = z_acc - z_h
            
            # The squared norm of the H-Acc vector.
            r  = dx*dx + dy*dy + dz*dz                        
            
            # The H-Don vector components in real coordinates.
            dx1 = x_don - x_h
            dy1 = y_don - y_h
            dz1 = z_don - z_h
            
            # The squared norm of the H-Don vector.
            r1 = dx1*dx1 + dy1*dy1 + dz1*dz1
                       
            # The (H_Acc,H-Don) cosinus is computed.
            cosa = (dx*dx1) + (dy*dy1) + (dz*dz1)
            cosa /=sqrt(r*r1)                        
            cosa = fmax(-1.0,fmin(1.,cosa))
            
            # The angle corresponding to the computed cosinus.
            angle = acos(cosa)
                        
            # If the angle is not in the H-Bond accepted interval, do not go futher.
            if angle < ang_min or angle > ang_max: continue
                                    
            ind = j + i*hbond_matrix.dimensions[1]
                        
            # The H bond matrix is updated with the new H bond found.
            hbond_mat[ind] = 1.0


