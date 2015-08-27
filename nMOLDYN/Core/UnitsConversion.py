"""This modules defines the conversion factors for the units used in nMOLDYN.
"""

# The units conversion factors.
# nMOLDYN standard units = THz (frequency), nm (distance), nm-1 (reciprocal distance), ps (time).
UNITS_CONV = {}

UNITS_CONV['frequency'] = (('thz'    , 1.0)               , \
                           ('cm-1'   , 33.356409519815209), \
                           ('mev'    , 4.1357319666954249), \
                           ('uev'    , 4135.7319666954245), \
                           ('rad s-1', 6.2831853071795862))

UNITS_CONV['distance'] = (('nm' , 1.0)   , \
                          ('ang', 1.0e1) , \
                          ('pm' , 1.0e3) , \
                          ('fm' , 1.0e6) , \
                          ('um' , 1.0e-3), \
                          ('mm' , 1.0e-6))

UNITS_CONV['q'] = (('nm-1' , 1.0)   , \
                   ('ang-1', 1.0e-1), \
                   ('pm-1' , 1.0e-3), \
                   ('fm-1' , 1.0e-6), \
                   ('um-1' , 1.0e3) , \
                   ('mm-1' , 1.0e6))

UNITS_CONV['time'] = (('ps', 1.0)   , \
                      ('fs', 1.0e3) , \
                      ('ns', 1.0e-3), \
                      ('us', 1.0e-6), \
                      ('ms', 1.0e-9))
