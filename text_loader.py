#!/usr/bin/env python3

def load_events(filename, energy_thr = 0.01, nof_events = -1):
    """
    load photon events from a text file (W,E,X,Y,Z,WX,WY,WZ), remove events with energy below threshold
    :param filename: file name to load events from
    :param energy_thr: energy threshold, events with energy below threshold will be thrown out
    """
    events = []

    k = 0
    with open(filename, "r+") as f:
        for line in f:

            line = line.strip()

            s = line.split()
            s = [x for x in s if x]  # remove empty lines

            if "GGG" not in s[0]:
                continue

            WT = float(s[1])
            E  = float(s[2])
            if E < energy_thr:
                break

            X = float(s[3])
            Y = float(s[4])
            Z = float(s[5])

            WX = float(s[6])
            WY = float(s[7])
            WZ = float(s[8])

            e = (W,E,X,Y,Z,WX,WY,WZ)
            events.append(e)

            k += 1
            if nof_events > 0:
                if k == nof_events:
                    break

    return events if len(events) > 0 else None
