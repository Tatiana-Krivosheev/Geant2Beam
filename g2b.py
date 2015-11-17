#!/usr/bin/env python

# all data are in mm, written out in cm

import sys
import math
import random
import struct
import numpy as np

def mm2cm(value):
    """
    Converts mm to cm
    """
    return value*0.1

def load_events(filename, energy_thr = 0.01):
    """
    load all events from a text file (W,E,X,Y,Z,WX,WY,WZ), remove events with energy below threshold
    :param filename: file name to load events from
    :param energy_thr: energy threshold, events with energy below threshold will be thrown out
    """
    events = []

    with open(filename) as f:
        for line in f:

            line = line.strip()
            s    = line.split()
            e    = []
            for n in s:
                e.append(float(n))

            if len(e) == 8:
                if e[1] >= energy_thr:
                    events.append(e)

    if len(events) == 0:
        return None

    return events

def shift_back(e, zshift):
    """
    Project back particles from current plane to plane shifted by zshift
    :param e: an event
    :param zshift: shift value, backward, in mm
    """

    if zshift == 0.0:
        return (e[2], e[3], e[4])

    x = e[2]
    y = e[3]
    z = e[4]
    
    wx = e[5]
    wy = e[6]
    wz = e[7]

    s = 0.0
    if wz != 0.0:
        s = zshift / wz
        
    return (x - wx*s, y - wy*s, z - wz*s)

    
def write_record_long(e, zshift, f, randomize = False):
    """
    Write MODE2 PhSF particle record

    :param e: an event
    :param zshift: shift value, backward, in mm
    :param f: file
    """

    LATCH = 8388608 # magic value, check PIRS-509 for description of bits in LATCH bitmask
    f.write(struct.pack("i", LATCH))

    E = -np.float32(e[1])
    f.write(struct.pack("f", E))

    if randomize:
        rotate_particle(e)
    
    (xx, yy, zz) = shift_back(e, zshift)

    X = np.float32(mm2cm(xx))
    f.write(struct.pack("f", X))

    Y = np.float32(mm2cm(yy))
    f.write(struct.pack("f", Y))

    U = np.float32(e[5])
    f.write(struct.pack("f", U))

    V = np.float32(e[6])
    f.write(struct.pack("f", V))

    WT = np.float32(e[0])
    if e[7] < 0.0: # keep W sign
        WT = np.float32(-WT)
    f.write(struct.pack("f", WT))

    ZLAST = np.float32(mm2cm(e[4] - 150.00)) # back by 150mm!
    f.write(struct.pack("f", ZLAST))

def rotate_2d( x, y, cs, sn ):
    """
    rotate x, y in he plane given angle sine and cosine
    """

    return (cs*x + sn*y, -sn*x + cs*y)

def random_angle():
    """
    sample random angle in azimuth, return tuple of (cos, sin)
    """
    phi = 2.0 * math.pi * random.random()
    cs  = math.cos( phi )
    sn  = math.sin( phi )

    return (cs, sn)

def rotate_particle(e):
    """
    Given particle phase space coordinates, rotate in XY plane
    """

    cs, sn = random_angle()

    x = e[2]
    y = e[3]

    e[2], e[3] = rotate_2d( x, y, cs, sn )
    
    wx = e[5]
    wy = e[6]

    e[5], e[6] = rotate_2d( wx, wy, cs, sn )    

def write_beam_long(header, events, zshift, filename, randomize):
    """
    Write BEAMnrc phase space file, MODE2 format
    """

    (mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP) = header

    with open(filename, "wb") as f:
        f.write(mode)
        f.write(struct.pack("i", NPPHSP))
        f.write(struct.pack("i", NPHOTPHSP))

        tmp = np.float32(EKMAX)
        f.write(struct.pack("f", tmp))

        tmp = np.float32(EKMIN)
        f.write(struct.pack("f", tmp))

        tmp = np.float32(NINCP)
        f.write(struct.pack("f", tmp))

        dummy = b"XXXXXXX" # 7 bytes to fill header up to 32bytes
        f.write(dummy)

        for e in events:
            write_record_long(e, zshift, f, randomize)

def make_header(nof_original, events):
    """
    Produce  from Geant data header for BeamNRC

    :param nof_original: number of original events used to get events in G4 phase space file
    :param events: all events from G4 phase space file
    """
    mode = b"MODE2"

    NPPHSP = len(events) # total nof of particle records

    NPHOTPHSP =  len(events) # total nof of photon records

    EKMAX = -1.0        # max kinetic energy, MeV
    EKMIN = 999999999.0 # min e- kinetic energy, MeV
    for e in events:
        energy = e[1]
        if energy > EKMAX:
            EKMAX = energy
        if energy < EKMIN:
            EKMIN = energy
            
    #EKMIN = 0.183944478631 # from previous PhSF
    
    NINCP = float(nof_original)

    return (mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP)

def main():
    """
    Process Geant phase space file and build BEAMnrc phase space file
    """

    if len(sys.argv) < 3:
        print("g2b input_fname output_phsf_name_without_extension <optional number of original decays, 10bil default>")
        return

    nof_decays = 10000000000
    try:
        n = int(sys.argv[3])
        nof_decays = n
    except:
        pass

    events = load_events(sys.argv[1])
    
    header = make_header(nof_decays, events)
    (mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP) = header

    print(mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP)

    write_beam_long(header, events, 0.0, sys.argv[2] + ".egsphsp1", False)
    
    return 0

if __name__ == '__main__':

    sys.exit(main())
