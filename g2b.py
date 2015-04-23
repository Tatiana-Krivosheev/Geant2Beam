#!/usr/bin/env python
#
# requires python 3.x. tested with 3.4.2
#

# all data are in mm, written out in cm

import math
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
                if e[1] > energy_thr:
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
    X = e[2]
    Y = e[3]
    Z = e[4]
    
    wx = e[5]
    wy = e[6]
    wz = e[7]
    
    s = 0.0
    if wz != 0.0:
        s = zshift / wz
        
    return (X-wx*s, Y-wy*s, Z-wz*s)

    
def write_record_long(e, zshift, f):
    """
    Write MODE2 PhSF particle record
    """

    LATCH = 8388608
    f.write(struct.pack("i", LATCH))

    E = -np.float32(e[1])
    f.write(struct.pack("f", E))
    
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


def write_beam_long(header, events, zshift, filename):
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
            write_record_long(e, zshift, f)
    
def make_header(nof_original, events):
    """
    Produce  from Geant data header for BeamNRC
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
            
    EKMIN = 0.183944478631 # from previous PhSF
    
    NINCP = float(nof_original)
    
    return (mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP)
    

events = load_events("../run25/photons")

header = make_header(10000000000, events)
(mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP) = header

print(mode, NPPHSP, NPHOTPHSP, EKMAX, EKMIN, NINCP)

write_beam_long(header, events, 9.98, "Q25.egsphsp1")

