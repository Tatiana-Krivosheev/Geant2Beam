#!/usr/bin/env python3

r"""
Convert from binary BEAM PhSF to text resembling Geant4 representation
"""

import sys
import math
import random
import struct
import numpy as np

import beam_loader

# all data are in cm, written out in mm
def cm2mm(value):
    """
    Converts cm to mm
    """
    return value*10.0


def write_g4_file(events, filename, Zsrc, Zdst):
    """
    Give list of events, print them in Geant4 format,
    move if necessary from Z source plane to Z destiantion plane
    """
    with open(filename, "wt+") as f:
        for event in events:
            wt = event[0]
            e  = event[1]
            x  = cm2mm(event[2])
            y  = cm2mm(event[3])
            zl = cm2mm(event[4])
            wx = event[5]
            wy = event[6]
            wz = event[7]

            (wt, e, x, y, z, wx, wy, wz) = move_event((wt, e, x, y, zl, wx, wy, wz), Zsrc, Zdst)

            l = "GGG:      {0}     {1}    {2}    {3}     {4}    {5}    {6}     {7}\n".format(wt,e,x,y,z,wx,wy,wz)
            f.write(l)


def move_event(event, Zsrc, Zdst):
    """
    Move event from Z source to Z destination
    """
    wt = event[0]
    e  = event[1]
    x  = event[2]
    y  = event[3]
    zl = event[4]
    wx = event[5]
    wy = event[6]
    wz = event[7]

    s = 0.0
    if wz != 0.0:
        s = (Zdst - Zsrc) / wz

    x = x + wx*s
    y = y + wy*s
    z = Zdst

    return (wt, e, x, y, z, wx, wy, wz)

def main(phsf, filename, Zsrc = 0.0, Zdst = 0.0):
    """
    read binary PhSF, output data into text format,
    shift from Z source to Z destination if necessary
    """
    events, nof_photons, nof_electrons, nof_positrons = beam_loader.load_events(phsf)
    print("{0} photons loaded, {1} electrons loaded, {2} positrons loaded".format(nof_photons, nof_electrons, nof_positrons))

    write_g4_file(events, filename)

    return 0

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("g2b input_phsp_fname output_G4_name")
        sys.exit(0)

    phsp_fname = sys.argv[1]

    g4_fname = sys.argv[2]

    rc = main(phsp_fname, g4_fname, 200.0, 197.5)

    sys.exit(rc)
