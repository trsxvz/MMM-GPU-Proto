#!/usr/bin/env python3
"""Tracks of a Perfetto trace: usage: pftrace-tracks.py <file.pftrace>

Lists what Perfetto (https://ui.perfetto.dev) will show for a trace written
by rocprofv3 (--output-format pftrace), without any library: the tracks
(processes, threads, GPU queues...), with their numbers of slices, their
categories and their most frequent names. The file is a sequence of
TracePacket messages of the protocol buffers of Perfetto; only the track
descriptors, the track events and the interned names are decoded."""
import sys
from collections import Counter, defaultdict


def varint(b, i):
    r = s = 0
    while True:
        c = b[i]
        i += 1
        r |= (c & 0x7F) << s
        if c < 0x80:
            return r, i
        s += 7


def fields(b):
    """(field number, wire type, value) of a message"""
    i, n = 0, len(b)
    while i < n:
        k, i = varint(b, i)
        f, w = k >> 3, k & 7
        if w == 0:
            v, i = varint(b, i)
        elif w == 1:
            v, i = b[i:i + 8], i + 8
        elif w == 2:
            m, i = varint(b, i)
            v, i = b[i:i + m], i + m
        elif w == 5:
            v, i = b[i:i + 4], i + 4
        else:
            raise ValueError(f"wire type {w}")
        yield f, w, v


def text(v):
    return bytes(v).decode("utf-8", "replace")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    with open(sys.argv[1], "rb") as f:
        data = memoryview(f.read())
    tracks = {}  # uuid -> [name, parent, kind]
    events = defaultdict(Counter)  # uuid -> names of the slices
    categories = defaultdict(Counter)
    kinds = Counter()
    names = defaultdict(dict)  # sequence -> iid -> name
    cats = defaultdict(dict)
    packets = 0
    for f, w, packet in fields(data):
        if f != 1 or w != 2:
            continue
        packets += 1
        sequence, descriptor, event, interned = 0, None, None, None
        for g, _, v in fields(packet):
            if g == 10:
                sequence = v
            elif g == 60:
                descriptor = v
            elif g == 11:
                event = v
            elif g == 12:
                interned = v
        if interned is not None:
            for g, w2, v in fields(interned):
                if g in (1, 2) and w2 == 2:
                    iid, name = 0, ""
                    for h, _, x in fields(v):
                        if h == 1:
                            iid = x
                        elif h == 2:
                            name = text(x)
                    (cats if g == 1 else names)[sequence][iid] = name
        if descriptor is not None:
            uuid, parent, name, kind = 0, None, "", "track"
            for g, w2, v in fields(descriptor):
                if g == 1:
                    uuid = v
                elif g == 5:
                    parent = v
                elif g == 2 and w2 == 2:
                    name = text(v)
                elif g in (3, 4) and w2 == 2:
                    kind = "process" if g == 3 else "thread"
                    pid = tid = None
                    label = ""
                    for h, w3, x in fields(v):
                        if h == 1 and w3 == 0:
                            pid = x
                        elif h == 2 and w3 == 0 and g == 4:
                            tid = x
                        elif w3 == 2 and ((g == 3 and h == 6) or
                                          (g == 4 and h == 5)):
                            label = text(x)
                    name = name or label
                    name += f" (pid {pid}" + (f", tid {tid})" if g == 4
                                              else ")")
                elif g == 8:
                    kind = "counter"
            # a descriptor may be written again without its name
            if name or uuid not in tracks:
                tracks[uuid] = [name, parent, kind]
        if event is not None:
            uuid, kind, name = 0, 0, None
            for g, w2, v in fields(event):
                if g == 11:
                    uuid = v
                elif g == 9:
                    kind = v
                elif g == 23 and w2 == 2:
                    name = text(v)
                elif g == 10 and w2 == 0:
                    name = names[sequence].get(v, f"name {v}")
                elif g == 22 and w2 == 2:
                    categories[uuid][text(v)] += 1
                elif g == 3 and w2 == 0:
                    categories[uuid][cats[sequence].get(v, f"category {v}")] += 1
            kinds[kind] += 1
            if kind in (1, 3) and name is not None:
                events[uuid][name] += 1
    print(f"{sys.argv[1]}: {len(data)} bytes, {packets} packets, "
          f"{len(tracks)} tracks, events: {kinds[1]} beginnings of slices, "
          f"{kinds[2]} ends, {kinds[3]} instants, {kinds[4]} counter values")

    def path(uuid):
        out = []
        while uuid in tracks and len(out) < 8:
            name, parent, kind = tracks[uuid]
            out.append(f"{name or '(unnamed)'} [{kind}]")
            uuid = parent
        return " < ".join(out)

    for uuid in sorted(set(tracks) | set(events),
                       key=lambda u: -sum(events[u].values())):
        n = sum(events[uuid].values())
        if n == 0 and uuid in tracks and tracks[uuid][2] != "counter":
            continue
        top = ", ".join(f"{k[:50]} x{c}" for k, c in
                        events[uuid].most_common(4))
        cat = ", ".join(f"{k} x{c}" for k, c in
                        categories[uuid].most_common(3))
        print(f"- {path(uuid) or f'track {uuid}'}: {n} slices"
              f"{' (' + cat + ')' if cat else ''}: {top}")


main()
