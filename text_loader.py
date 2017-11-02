#!/usr/bin/env python3

def load_events(filename, energy_thr = 0.01):
    """
    load all events from a text file (W,E,X,Y,Z,WX,WY,WZ), remove events with energy below threshold
    :param filename: file name to load events from
    :param energy_thr: energy threshold, events with energy below threshold will be thrown out
    """
    events = []

    with open(filename, "r+") as f:
        for line in f:

            line = line.strip()

            s = line.split()
            s = [x for x in s if x]  # remove empty lines

            WT = float(s[0])
            E  = float(s[1])
            if E < energy_thr:
                break

            X = float(s[2])
            Y = float(s[3])
            Z = float(s[4])

            WX = float(s[5])
            WY = float(s[6])
            WZ = float(s[7])

            e = (W,E,X,Y,Z,WX,WY,WZ)
            events.append(e)

    return events if len(events) > 0 else None
