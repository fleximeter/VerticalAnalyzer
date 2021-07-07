"""
File: vslice.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains the v_slice class for vertical slicing with music21.
Copyright (c) 2021 by Jeff Martin.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import music21
from statistics import pstdev


class VSlice:
    def __init__(self, duration=1, quarter_duration=1, measure=None):
        """
        Creates a v_slice
        :param duration: The duration of the slice, in seconds
        """
        self._cardinality = 0  # The cardinality of the slice
        self._chord = None  # The chord
        self._duration = duration  # The duration of the slice in seconds
        self._forte = None  # The Forte name of the chord
        self._ins = None  # The INS of the slice
        self._ipseg = []  # The ipseg of the slice
        self._lns = None  # The LNS of the slice
        self._lower_bound = None  # The lower bound of the slice.
        self._measure = measure  # The measure number in which the slice begins
        self._mt = None  # The MT of the slice
        self._meant = None  # The mT of the slice
        self._ns = None  # The NS of the slice. If the lower and upper bounds are defined, but LNS and UNS are not,
        # the NS represents the entire pitch area encompassed by the piece. Otherwise it is None.
        self._pcset = []  # The pcset of the current v_slice
        self._pitch_list = []  # A list of pitches in the chord
        self._pitches = set()  # A list of distinct p-integers represented in the slice
        self._pitches_sorted = []  # A sorted list of pitches
        self._prime_form = None  # The prime form of the slice
        self._ps = None  # The PS of the slice
        self._quarter_duration = quarter_duration  # The duration in quarters
        self._uns = None  # The UNS of the slice
        self._upper_bound = None  # The upper bound of the slice.

    @property
    def cardinality(self):
        return self._cardinality

    @property
    def chord(self):
        return self._chord

    @property
    def duration(self):
        return self._duration

    @duration.setter
    def duration(self, value):
        self._duration = value

    @property
    def forte(self):
        return self._forte

    @property
    def ins(self):
        return self._ins

    @property
    def ipseg(self):
        return self._ipseg

    @property
    def lns(self):
        return self._lns

    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, value):
        self._lower_bound = value

    @property
    def measure(self):
        return self._measure

    @property
    def meant(self):
        return self._meant

    @property
    def mt(self):
        return self._mt

    @property
    def ns(self):
        return self._ns

    @property
    def pcset(self):
        return self._pcset

    @property
    def pitch_list(self):
        return self._pitch_list

    @property
    def pitches(self):
        return self._pitches

    @property
    def pitches_sorted(self):
        return self._pitches_sorted

    @property
    def prime_form(self):
        return self._prime_form

    @property
    def ps(self):
        return self._ps

    @property
    def quarter_duration(self):
        return self._quarter_duration

    @quarter_duration.setter
    def quarter_duration(self, value):
        self._quarter_duration = value

    @property
    def uns(self):
        return self._uns

    @property
    def upper_bound(self):
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, value):
        self._upper_bound = value

    def add_pitches(self, pitches, pitch_names=None):
        """
        Adds pitches to the v_slice
        :param pitches: A collection of pitches to add
        """
        # Add each pitch to the chord
        for p in pitch_names:
            self._pitch_list.append(p)
        self._cardinality = len(self._pitch_list)

        for pitch in pitches:
            self._pitches.add(pitch)

    def calculate_meant(self):
        """
        Calculates the mean trajectory
        """
        if self._ps > 0 and self._upper_bound != None and self._lower_bound != None:
            mean = 0
            for p in self._pitches:
                mean += p
            mean /= self._ps
            self._meant = mean - (self._upper_bound + self._lower_bound) / 2

    def get_ipseg_stdev(self):
        """
        Gets the population standard deviation of the v_slice ipseg
        :return: The standard deviation
        """
        if len(self._ipseg) > 0:
            return pstdev(self._ipseg)
        else:
            return None

    def get_ipseg_string(self):
        """
        Gets the ipseg as a string
        :return: The ipseg as a string
        """
        ipseg = "\"<"
        for ip in self._ipseg:
            ipseg += str(ip) + ", "
        if ipseg[len(ipseg) - 1] == " ":
            ipseg = ipseg[:-2]
        ipseg += ">\""
        return ipseg

    def run_calculations(self):
        """
        Calculates information about the v_slice. You must set the lower and upper bounds before running this
        method. You should also combine any v_slices that you want to combine before running this method,
        to avoid making unnecessary computations.
        :return: None
        """
        # Make a sorted pitch list
        self._pitches_sorted = list(self._pitches)
        self._pitches_sorted.sort()
        for i in range(1, len(self._pitches_sorted)):
            self._ipseg.append(self._pitches_sorted[i] - self._pitches_sorted[i - 1])

        # Calculate ps and ins
        self._ps = len(self._pitches)
        if self._ps > 0:
            self._ins = self._pitches_sorted[len(self._pitches_sorted) - 1] - self._pitches_sorted[0] + 1 - self._ps
        else:
            self._ins = 0

        # Calculate uns, lns, ns, and mt
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._ps is None or self._ps == 0:
                self._lns = None
                self._uns = None
                self._mt = None
                self._ns = self._upper_bound - self._lower_bound + 1
            else:
                self._lns = self._pitches_sorted[0] - self._lower_bound
                self._uns = self._upper_bound - self._pitches_sorted[len(self._pitches_sorted) - 1]
                self._mt = (self._lns - self._uns) / 2

        # Create music21 Chord object to represent the v_slice
        self._chord = music21.chord.Chord()
        for p in self._pitch_list:
            self._chord.add(p)

        # Calculate set theory info
        self._forte = self._chord.forteClass
        self._prime_form = "["
        self._pcset = "{"
        pcset = []
        for pc in self._chord.primeForm:
            if pc == 10:
                self._prime_form += "A"
            elif pc == 11:
                self._prime_form += "B"
            else:
                self._prime_form += str(pc)
        self._prime_form += "]"
        for p in self._pitches:
            pc = p % 12
            if pc < 0:
                pc += 12
            pcset.append(pc)
        pcset.sort()
        for pc in pcset:
            if pc == 10:
                self._pcset += "A"
            elif pc == 11:
                self._pcset += "B"
            else:
                self._pcset += str(pc)
        self._pcset += "}"
