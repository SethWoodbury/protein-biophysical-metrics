# SPDX-License-Identifier: BSD-2-Clause
#
# AGGRESCAN numerical conventions adapted from psico.aggrescanning.
# Copyright (c) 2022, Thomas Holder
# Copyright (c) 2022, Leukocare AG
#
# BSD 2-Clause License
#
# Copyright (c) 2010-2019, Thomas Holder
# Copyright (c) 2010-2012, Steffen Schmidt
# Copyright (c) 2005-2009, Ezequiel Panepucci (see psico/helping.py)
# Copyright (c) 2011, Tsjerk Wassenaar (see psico/querying.py)
# Copyright (c) 2009, Sean Law (see psico/exporting.py)
# Copyright (c) 2011, Takanori Nakane (spectrum_states)
# Copyright (c) 2006-2012, Christoph Gohlke (snippet in angle_between_domains)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Sequence aggregation descriptors using psico's AGGRESCAN conventions.

Reference implementation and BSD-2-Clause notice from pymol-psico commit
f42cb1e4f565685b95fba63820b5f37db32b5954, psico/aggrescanning.py and LICENSE.
Method DOI 10.1186/1471-2105-8-65.

The historical web-compatible profile uses mean anionic and cationic
scores for the N- and C-terminal padding, respectively. Off-center positions
reuse the nearest window average. A hotspot has at least five consecutive
non-Pro residues above the threshold. Its area sums the profile excess over
the threshold, with an extra half-height term when the run starts at the
first residue. The smallest padded window needs three actual residues.
"""

from collections.abc import Iterator

AGGRESCAN_GROUP = "aggrescan"

# Sánchez de Groot et al. (2005), DOI 10.1186/1472-6807-5-18.
_A3V = {
    "I": 1.822,
    "F": 1.754,
    "V": 1.594,
    "L": 1.380,
    "Y": 1.159,
    "W": 1.037,
    "M": 0.910,
    "C": 0.604,
    "A": -0.036,
    "T": -0.159,
    "S": -0.294,
    "P": -0.334,
    "G": -0.535,
    "K": -0.931,
    "H": -1.033,
    "Q": -1.231,
    "R": -1.240,
    "N": -1.302,
    "E": -1.412,
    "D": -1.836,
}
_N_TERMINAL_A3V = (_A3V["D"] + _A3V["E"]) / 2
_C_TERMINAL_A3V = (_A3V["K"] + _A3V["R"]) / 2
_HOTSPOT_THRESHOLD = -0.02
_MIN_HOTSPOT_LENGTH = 5
_WINDOWS = ((75, 5), (175, 7), (300, 9))


def _hotspots(sequence: str, profile: list[float]) -> Iterator[tuple[int, int, float]]:
    start = 0
    area = 0.0
    # The extra Pro closes a hotspot reaching the end of the sequence.
    for end, residue in enumerate(sequence + "P"):
        if residue != "P" and (height := profile[end] - _HOTSPOT_THRESHOLD) > 0:
            area += height
            continue
        if end - start >= _MIN_HOTSPOT_LENGTH:
            if start == 0:
                area += (profile[0] - _HOTSPOT_THRESHOLD) / 2
            yield start + 1, end, area
        start = end + 1
        area = 0.0


def describe_aggrescan(
    sequence: str, *, details: dict | None = None,
) -> dict[str, float | None]:
    """Score an uppercase canonical sequence and optionally record its hotspots.

    Na4vSS is 100 times the average a4v profile and THSAr is total hotspot
    area per residue. A one- or two-residue sequence has no valid window and
    both outputs are None. Reported hotspot intervals are 1-based inclusive.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    if not sequence:
        raise ValueError("sequence is empty")
    try:
        a3v = [_A3V[residue] for residue in sequence]
    except KeyError as error:
        raise ValueError(f"sequence contains noncanonical residue code {error.args[0]!r}") from error

    length = len(sequence)
    window = next((width for limit, width in _WINDOWS if length <= limit), 11)
    if length + 2 < window:
        if details is not None:
            details["aggrescan"] = {
                "window_size": window,
                "hotspots": None,
                "reason": "The padded window requires at least three residues.",
            }
        return {"aggrescan_na4vss": None, "aggrescan_thsar": None}

    padded = [_N_TERMINAL_A3V, *a3v, _C_TERMINAL_A3V]
    averages = [
        sum(padded[start:start + window]) / window
        for start in range(len(padded) - window + 1)
    ]
    off_center = window // 2 - 1
    profile = [averages[0]] * off_center + averages + [averages[-1]] * off_center

    intervals = [] if details is not None else None
    total_area = 0.0
    for start, end, area in _hotspots(sequence, profile):
        total_area += area
        if intervals is not None:
            intervals.append({"start": start, "end": end, "area": area})
    if details is not None:
        details["aggrescan"] = {"window_size": window, "hotspots": intervals}
    return {
        "aggrescan_na4vss": sum(profile) / length * 100,
        "aggrescan_thsar": total_area / length,
    }
