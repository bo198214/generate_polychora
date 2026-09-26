#!/usr/bin/env python3
"""
Example sounds from the modal analysis of the polychora (modal_output/*.json).

Every polychoron is struck once and rings as a sum of exponentially decaying sinusoids,
one per mode multiplet g (modal synthesis):

    v(t) = sum_g gain_g * exp(-pi * loss * f_g * t) * cos(2 pi f_g t)
    y    = v convolved with the mallet's force pulse (half sine of the contact time)

v is the velocity at the hit point after a unit impulse (--pickup changes that), gain_g
the multiplet's excitation for the chosen kind of hit (modal_output "gains"), and a
constant loss factor makes every mode ring for the same number of cycles
(T60 = 2.2 / (loss * f)).

Common size scaling: all polychora get the same circumradius R. The modal files hold
dimensionless frequencies for R = 1 and shear wave speed c_s = 1, so one factor for all,

    f_Hz = scale * frequency,   scale = c_s / R,

puts them into the audible range. (For one shape every size measure is proportional;
across shapes the circumradius predicts the fundamental best: ~x1.3 spread at equal R,
~x2.2 at equal volume, x22 at equal edge length.)

Symmetry breaking (--detune): a real object is never perfectly symmetric, so each
multiplet of multiplicity m splits into m partials with slightly different frequencies
(the beating of bells). Random-matrix model of an unknown small imperfection: the relative
detunings are the eigenvalues of a random symmetric m x m matrix (GOE, level repulsion),
scaled to the given RMS; the multiplet's gain is divided among the partials with
Porter-Thomas weights z_k^2 / |z|^2, z ~ N(0, I), so the attack is unchanged. The
frequencies are fixed per object (seed from --seed and the name), the weights per hit.

Usage:
  python synth_modal.py                        # all -> modal_output/sounds/*.wav + _tour.wav
  python synth_modal.py --scale 2000 pen ex    # other scaling, selected polychora
  python synth_modal.py --site vertex          # strike a vertex instead of a typical point
  python synth_modal.py --site all             # typical hit, each cell class, vertex
  python synth_modal.py --detune 0.003 --out-dir modal_output/sounds_detune
"""
import argparse
import glob
import json
import os
import sys
import zlib

import numpy as np
from scipy.io import wavfile
from scipy.signal import fftconvolve


def strikes(doc, site):
    """(label, gains) of the hits to render: mean = uniformly random point on the boundary."""
    g = doc["gains"]
    if site == "mean":
        return [("typical hit", g["mean"])]
    if site == "vertex":
        return [("vertex", g["vertex"])]
    if site == "cell":
        s = max(g["cells"], key=lambda s: s["count"])
        return [(s["label"], s)]
    return ([("typical hit", g["mean"])] + [(s["label"], s) for s in g["cells"]]
            + [("vertex", g["vertex"])])


def frequencies_hz(doc, scale):
    """Multiplet frequencies in Hz for scale = c_s / R."""
    return np.asarray(doc["modes"]["frequency"]) * scale


def partials(doc, gain, args, hit=0):
    """Frequencies [Hz] and gains of the partials of one hit. Without --detune one per
    multiplet; with it, m per multiplet: the object's imperfection fixes the split
    frequencies, the hit (index `hit`) draws how the multiplet's gain divides among them."""
    freq, gain = frequencies_hz(doc, args.scale), np.asarray(gain)
    if args.detune <= 0:
        return freq, gain
    key = zlib.crc32(doc["name"].encode())
    obj = np.random.default_rng([args.seed, key])
    rng = np.random.default_rng([args.seed, key, hit + 1])
    fs, gs = [], []
    for f, m, g in zip(freq, doc["modes"]["multiplicity"], gain):
        A = obj.standard_normal((m, m))
        lam = np.linalg.eigvalsh((A + A.T) / 2) / np.sqrt((m + 1) / 2)   # unit RMS
        z = rng.standard_normal(m)
        fs.append(f * (1 + args.detune * lam))
        gs.append(g * z * z / (z @ z))
    return np.concatenate(fs), np.concatenate(gs)


def ring(freq, gain, args, seconds):
    """One strike: modal response at the hit point, filtered by the mallet pulse."""
    use = freq < min(20000.0, 0.45 * args.sr)
    t = np.arange(int(seconds * args.sr)) / args.sr
    y = np.zeros_like(t)
    if not use.any():
        return y, 0.0
    for f, g in zip(freq[use], gain[use]):
        w = 2 * np.pi * f
        env = g * np.exp(-np.pi * args.loss * f * t)
        if args.pickup == "displacement":
            y += env / w * np.sin(w * t)
        elif args.pickup == "velocity":
            y += env * np.cos(w * t)
        else:
            y -= env * w * np.sin(w * t)
    # half-sine force pulse; auto: its first spectral zero (at 1.5 / contact) sits on the
    # highest exported mode, so the truncation of the modal set is not audible
    contact = args.contact_ms / 1000 if args.contact_ms else 1.5 / freq[use].max()
    n = max(1, round(contact * args.sr))
    pulse = np.sin(np.pi * (np.arange(n) + 0.5) / n)
    return fftconvolve(y, pulse / pulse.sum())[:len(t)], contact


def finish(y, fade, sr):
    """Raised-cosine fade-out and peak normalisation to -1 dBFS."""
    n = min(len(y), int(fade * sr))
    y[len(y) - n:] *= 0.5 * (1 + np.cos(np.linspace(0, np.pi, n)))
    peak = np.abs(y).max()
    return y * (0.891 / peak) if peak > 0 else y


def write_wav(path, y, sr):
    wavfile.write(path, sr, np.round(y * 32767).astype(np.int16))


def main():
    ap = argparse.ArgumentParser(description="Example sounds from modal_output/*.json")
    ap.add_argument("names", nargs="*", help="polychora (default: all in --modal-dir)")
    ap.add_argument("--scale", type=float, default=1000.0,
                    help="common size scaling c_s/R in Hz, all polychora with the same "
                         "circumradius R (default 1000, e.g. steel with R = 3.2 m)")
    ap.add_argument("--site", choices=("mean", "cell", "vertex", "all"), default="mean",
                    help="mean: typical hit (random point); cell: most frequent cell class; "
                         "vertex; all: mean, every cell class, vertex in sequence")
    ap.add_argument("--direction", choices=("normal", "tangential"), default="normal")
    ap.add_argument("--pickup", choices=("displacement", "velocity", "acceleration"),
                    default="velocity", help="signal = this quantity at the hit point")
    ap.add_argument("--loss", type=float, default=0.002,
                    help="material loss factor (steel ~1e-4..1e-3, glass ~1e-3, wood ~1e-2)")
    ap.add_argument("--contact-ms", type=float, default=0.0,
                    help="mallet contact time in ms (default: adapted to each polychoron's modes)")
    ap.add_argument("--detune", type=float, default=0.0,
                    help="RMS relative detuning of degenerate modes by imperfections (0 = ideal "
                         "object; 0.0002..0.002 bell-like, 0.003 clear slow beating, 0.01 chorus)")
    ap.add_argument("--seed", type=int, default=0, help="which imperfect exemplar (with --detune)")
    ap.add_argument("--duration", type=float, default=4.0, help="seconds per strike")
    ap.add_argument("--tour-seconds", type=float, default=2.5,
                    help="seconds per polychoron in _tour.wav (0 = no tour)")
    ap.add_argument("--sr", type=int, default=48000)
    ap.add_argument("--modal-dir", default="modal_output")
    ap.add_argument("--out-dir", default=os.path.join("modal_output", "sounds"))
    args = ap.parse_args()

    docs = {}
    for p in sorted(glob.glob(os.path.join(args.modal_dir, "*.json"))):
        with open(p, encoding="utf-8") as f:
            doc = json.load(f)
        docs[doc["name"]] = doc
    names = args.names or sorted(docs)
    missing = [n for n in names if n not in docs]
    if missing:
        sys.exit(f"no modal data for {missing} in {args.modal_dir}/ (run modal_analysis.py)")
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"scale {args.scale:g} Hz (c_s/R), loss {args.loss:g}, site {args.site}, "
          f"pickup {args.pickup}, detune {args.detune:g} -> {args.out_dir}/")
    rows = []
    for name in names:
        doc = docs[name]
        freq = frequencies_hz(doc, args.scale)
        parts, labels = [], []
        for k, (label, g) in enumerate(strikes(doc, args.site)):
            y, contact = ring(*partials(doc, g[args.direction], args, k), args, args.duration)
            parts.append(finish(y, 0.5, args.sr))
            labels.append(label)
        write_wav(os.path.join(args.out_dir, name + ".wav"), np.concatenate(parts), args.sr)
        rows.append((name, freq[0], freq[-1]))
        print(f"  {name:10s} f1 = {freq[0]:7.1f} Hz  top mode {freq[-1]:7.1f} Hz  "
              f"{len(freq):3d} multiplets  T60(f1) = {2.2 / (args.loss * freq[0]):5.1f} s  "
              f"contact {contact * 1000:.2f} ms  [{', '.join(labels)}]")

    if args.tour_seconds > 0 and len(names) > 1:
        legend, segs = [], []
        site = "mean" if args.site == "all" else args.site
        for i, (name, f1, _) in enumerate(sorted(rows, key=lambda r: -r[1])):
            doc = docs[name]
            y, _ = ring(*partials(doc, strikes(doc, site)[0][1][args.direction], args), args,
                        args.tour_seconds)
            segs.append(finish(y, 0.4, args.sr))
            legend.append(f"{i * args.tour_seconds:6.1f} s  {name:10s} {f1:7.1f} Hz  "
                          f"{doc['description']}")
        write_wav(os.path.join(args.out_dir, "_tour.wav"), np.concatenate(segs), args.sr)
        with open(os.path.join(args.out_dir, "_tour.txt"), "w", encoding="utf-8") as f:
            f.write(f"scale {args.scale:g} Hz (c_s/R, same circumradius), strike: {site}, "
                    f"detune {args.detune:g}, {args.tour_seconds:g} s each, "
                    f"highest to lowest fundamental\n")
            f.write("\n".join(legend) + "\n")
        print(f"  _tour.wav: {len(segs)} polychora, highest to lowest (legend in _tour.txt)")

    lo = min(r[1] for r in rows)
    hi = max(r[2] for r in rows)
    print(f"fundamentals {lo:.1f} .. {max(r[1] for r in rows):.1f} Hz, all modes up to {hi:.1f} Hz")
    if lo < 20 or hi > 20000:
        print("WARNING: part of the spectrum lies outside the audible range 20 Hz .. 20 kHz; "
              "adjust --scale")


if __name__ == "__main__":
    main()
