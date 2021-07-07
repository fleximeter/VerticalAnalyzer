"""
File: analyze.py
Author: Jeff Martin
Email: jeffreymartin@outlook.com
This file contains functions for analyzing vertical slices.
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
from vslice import VSlice
from fractions import Fraction


def analyze(section_name, input_xml, output_path, general_path="", general_command="", first=-1, last=-1, uselocal="N",
            output_measure=False):
    """
    Performs a vertical analysis on the given stream and writes a report to CSV
    :param input_xml: The musicxml file to analyze
    :param output_path: The file path for the report
    """
    stream = music21.converter.parse(input_xml)
    parts = []
    for item in stream:
        if type(item) == music21.stream.Part:
            parts.append(item)
    slices = slice_parts(parts, get_slice_num(parts), first, last, uselocal, output_measure)
    write_report(output_path, slices, general_path, general_command, section_name)


def clean_slices(slices):
    """
    Cleans up a list of v_slices
    :param slices: A list of v_slices
    """
    # Remove duplicate slices, and update durations
    if len(slices) > 0:
        i = 1
        while i < len(slices):
            if slice_compare(slices[i - 1], slices[i]):
                slices[i - 1].duration += slices[i].duration
                slices[i - 1].quarter_duration += slices[i].quarter_duration
                del slices[i]
            else:
                i += 1


def factor(n):
    """
    Factors a positive integer
    :param n: An integer
    :returns: A list of factors, in sorted order, including duplicates
    """
    factors = []
    d = 1
    while d <= int(n ** 0.5):
        if n % d == 0:
            factors.append(d)
            n //= d
        else:
            d += 1
        if d == 1:
            d += 1
        if d > int(n ** 0.5):
            factors.append(n)
    factors.sort()
    return factors


def get_bounds(slices):
    """
    Gets the upper and lower bounds of a list of v_slices
    :param slices: A list of v_slices
    :return: The lower and upper bounds as a tuple. The lower bound is index 0,
    and the upper bound is index 1. If the slices contain no pitches, each of
    the bounds will be None.
    """

    lower_bound = None
    upper_bound = None

    for i in range(0, len(slices)):
        if (lower_bound is None or upper_bound is None) and len(slices[i].pitches_sorted) > 0:
            lower_bound = slices[i].pitches_sorted[0]
            upper_bound = slices[i].pitches_sorted[len(slices[i].pitches_sorted) - 1]
        if len(slices[i].pitches_sorted) > 0 and lower_bound is not None and upper_bound is not None:
            if slices[i].pitches_sorted[0] < lower_bound:
                lower_bound = slices[i].pitches_sorted[0]
            if slices[i].pitches_sorted[len(slices[i].pitches_sorted) - 1] > upper_bound:
                upper_bound = slices[i].pitches_sorted[len(slices[i].pitches_sorted) - 1]

    return lower_bound, upper_bound


def get_piece_bounds(parts):
    """
    Determines the lower and upper bounds of a piece
    :param parts: A list of parts
    :return: The lower and upper bounds as a tuple.
    """
    lower = None
    upper = None
    for part in parts:
        for item in part:
            if type(item) == music21.stream.Measure:
                for item2 in item:
                    if type(item2) == music21.stream.Voice:
                        for item3 in item2:
                            if type(item3) == music21.note.Note or type(item3) == music21.chord.Chord:
                                for pitch in item3.pitches:
                                    if lower is None:
                                        lower = pitch.midi - 60
                                    if upper is None:
                                        upper = pitch.midi - 60
                                    if lower > pitch.midi - 60:
                                        lower = pitch.midi - 60
                                    if upper < pitch.midi - 60:
                                        upper = pitch.midi - 60
                    elif type(item2) == music21.note.Note or type(item2) == music21.chord.Chord:
                        for pitch in item2.pitches:
                            if lower is None:
                                lower = pitch.midi - 60
                            if upper is None:
                                upper = pitch.midi - 60
                            if lower > pitch.midi - 60:
                                lower = pitch.midi - 60
                            if upper < pitch.midi - 60:
                                upper = pitch.midi - 60

    return lower, upper


def get_slice_num(parts):
    """
    Determines the number of slices per quarter note based on subdivisions of the note.
    :param parts: A stream of parts
    :returns: The number of slices per quarter note
    """
    # A collection of all the unique denominators we find
    denominators = {}
    denominators_list = []

    # Find all the unique denominators
    for part in parts:
        for stream in part:
            if type(stream) == music21.stream.Measure:
                for item in stream:
                    if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(item) == \
                            music21.chord.Chord:
                        ql = item.duration.quarterLength
                        if type(item.duration.quarterLength) != Fraction:
                            ql = Fraction(item.duration.quarterLength)
                        if ql.denominator not in denominators:
                            denominators[ql.denominator] = True

    # Get the LCM and return it. This is the number of slices per quarter note that we need.
    for item in denominators:
        denominators_list.append(item)
    return lcm(denominators_list)


def lcm(integers):
    """
    Computes the LCM of a list of positive integers
    :param integers: A list of positive integers
    :return: The LCM
    """
    factors = {}  # A dictionary of individual factors and their multiplicities
    multiple = 1  # The LCM

    for num in integers:
        cur_factors = factor(num)  # The factors of the current number
        current = 1  # The current factor we are considering
        count = 0  # The number of occurrences of that factor
        for i in range(len(cur_factors)):
            # If we found another occurrence of that factor, increase the count
            if cur_factors[i] == current:
                count += 1
            # Otherwise record the count and move on
            else:
                if current not in factors:
                    factors[current] = count
                elif factors[current] < count:
                    factors[current] = count
                current = cur_factors[i]
                count = 1
            # If we are done, record the count of the last factor
            if i + 1 == len(cur_factors):
                if current not in factors:
                    factors[current] = count
                elif factors[current] < count:
                    factors[current] = count

    # Compute the LCM
    for item in factors:
        multiple *= item ** factors[item]
    # print(multiple)
    return multiple


def set_slice_bounds(slices, bounds):
    """
    Sets the bounds of a list of v_slices
    :param slices: A list of v_slices
    """

    for i in range(len(slices)):
        slices[i].lower_bound = bounds[0]
        slices[i].upper_bound = bounds[1]


def slice_parts(parts, n=1680, first=-1, last=-1, uselocal="N", output_numbers=False):
    """
    Takes n vertical slices of each beat from each of the parts. Note that beats are always quarter notes
    in music21. The parts do not need to have the same time signature for each measure: each slice is taken
    independently of the time signature. The parts do not even need to have the same number of total beats.
    However, it is assumed that a quarter note in any given part is equal in duration to a quarter note in
    any other part (this means that all parts must share the same tempo for a quarter note).
    :param parts: A list of parts
    :return: A list of v_slices
    """

    final_slices = []  # Holds the finalized slices to return
    n = Fraction(n, 1)  # The number of slices per quarter note
    transpose = [0 for i in range(len(parts))]  # The amount by which to transpose, for each part
    next_indices = [0 for i in range(len(parts))]  # The index of the next measure, for each part
    next_measure = -1  # The number of the next measure
    tempo = 60.0  # We assume a tempo of 60 to begin
    tempo_multiplier = 10  # This is in place to avoid floats

    if len(parts) == 0:
        print("No parts were provided")

    else:
        # Determine the index of the first measure in each part. a is the part index,
        # and b is the index of the item inside the current part (which may or may not be a measure)
        for a in range(len(parts)):
            found_first = False
            for b in range(len(parts[a])):
                if type(parts[a][b]) == music21.stream.Measure:
                    if parts[a][b].number >= first:
                        next_measure = parts[a][b].number
                        next_indices[a] = b
                        found_first = True
                # No need to continue after the first measure was found
                if found_first:
                    break

        # We consider each measure separately. When we have finished the last measure,
        # the next_measure will reset to -1 and we will stop.
        while next_measure != -1:
            if output_numbers:
                print("Measure", next_measure)
            # The slices taken for this measure
            measure_slices = []

            # Consider each part separately for this measure
            for a in range(len(parts)):
                # Tracks the number of slices taken for the current part in the current measure
                num_slices_taken = 0
                for item in parts[a][next_indices[a]]:
                    last_item_was_voice = False
                    furthest_voice_slice = 0

                    # MusicXML doesn't handle transposition properly for 8va and 8vb clefs, so we need manual
                    # transposition. Record for the future.
                    if type(item) == music21.clef.Bass8vaClef or type(item) == music21.clef.Treble8vaClef:
                        transpose[a] = 12
                    elif type(item) == music21.clef.Bass8vbClef or type(item) == music21.clef.Treble8vbClef:
                        transpose[a] = -12
                    elif isinstance(item, music21.clef.Clef):
                        transpose[a] = 0

                    # Update the tempo if we find a new one
                    if type(item) == music21.tempo.MetronomeMark:
                        tempo = item.number
                        tempo_multiplier = 10 ** str(tempo)[::-1].find(".")

                    # If we have found multiple voices in the same part in the same measure
                    if type(item) == music21.stream.Voice:
                        last_item_was_voice = True

                        # Track the start point for the voice
                        slice_start = num_slices_taken

                        for item2 in item:
                            # We can only take slices of notes, rests, or chords
                            if type(item2) == music21.note.Note or type(item2) == music21.note.Rest or type(
                                    item2) == music21.chord.Chord:
                                ql = item2.duration.quarterLength
                                if type(item2.duration.quarterLength) != Fraction:
                                    ql = Fraction(item2.duration.quarterLength)

                                num_slices = int(ql * n)
                                # the pitches are considered as integers in p-space. The p_names hold pitch names
                                # which is often more convenient for humans.
                                pitches_in_item = []
                                p_names_in_item = []

                                # We use Morris's p-space. Obviously rests do not have pitches.
                                if type(item2) != music21.note.Rest:
                                    for p in item2.pitches:
                                        name = p.name
                                        octave = p.octave + (transpose[a] // 12)
                                        pitches_in_item.append(p.midi - 60 + transpose[a])
                                        p_names_in_item.append(name + str(octave))

                                # Perform slicing. num_slices is the number of slices we take for the current object.
                                for j in range(num_slices):
                                    if num_slices_taken >= len(measure_slices):
                                        measure_slices.append(
                                            VSlice(Fraction(60 * tempo_multiplier, int(tempo * tempo_multiplier) * n),
                                                   Fraction(1, n.numerator), parts[a][next_indices[a]].number))
                                    measure_slices[num_slices_taken].add_pitches(pitches_in_item, p_names_in_item)
                                    num_slices_taken += 1

                        # Record the furthest slice reached in this voice if necessary
                        if furthest_voice_slice < num_slices_taken:
                            furthest_voice_slice = num_slices_taken

                        # Reset the slice counter to start on the next voice
                        num_slices_taken = slice_start

                    # If we just evaluated a voice and are done, we need to reset the slice counter
                    elif last_item_was_voice and num_slices_taken < furthest_voice_slice:
                        num_slices_taken = furthest_voice_slice

                    # We can only take slices of notes, rests, or chords
                    if type(item) == music21.note.Note or type(item) == music21.note.Rest or type(
                            item) == music21.chord.Chord:
                        ql = item.duration.quarterLength
                        if type(item.duration.quarterLength) != Fraction:
                            ql = Fraction(item.duration.quarterLength)

                        num_slices = int(ql * n)

                        # the pitches are considered as integers in p-space. The p_objects hold pitch names which
                        # is often more convenient for humans.
                        pitches_in_item = []
                        p_names_in_item = []

                        # We use Morris's p-space. Obviously rests do not have pitches.
                        if type(item) != music21.note.Rest:
                            for p in item.pitches:
                                name = p.name
                                octave = p.octave + (transpose[a] // 12)
                                pitches_in_item.append(p.midi - 60 + transpose[a])
                                p_names_in_item.append(name + str(octave))

                        # Perform slicing. num_slices is the number of slices we take for the current object.
                        for j in range(num_slices):
                            if num_slices_taken >= len(measure_slices):
                                measure_slices.append(
                                    VSlice(Fraction(60 * tempo_multiplier, int(tempo * tempo_multiplier) * n),
                                           Fraction(1, n.numerator), parts[a][next_indices[a]].number))
                            measure_slices[num_slices_taken].add_pitches(pitches_in_item, p_names_in_item)
                            num_slices_taken += 1

            # Clean up the slices from this measure
            clean_slices(measure_slices)
            for item in measure_slices:
                final_slices.append(item)

            # Find the next measure for each part
            for a in range(len(parts)):  # a is the part index
                found_next = False
                next_measure = -1
                # We start at the item after the current measure
                for b in range(next_indices[a] + 1, len(parts[a])):
                    if type(parts[a][b]) == music21.stream.Measure:
                        next_measure = parts[a][b].number
                        next_indices[a] = b
                        found_next = True
                    # No need to continue after the first measure was found
                    if found_next:
                        break

            # If we've analyzed the last measure, it's time to stop analyzing
            if next_measure > last > -1:
                next_measure = -1

    # Cleanup and return
    clean_slices(final_slices)
    bounds = None
    if uselocal == "N":
        bounds = get_piece_bounds(parts)
    else:
        bounds = get_bounds(final_slices)
    set_slice_bounds(final_slices, bounds)
    for f_slice in final_slices:
        f_slice.run_calculations()
        f_slice.calculate_meant()
    return final_slices


def slice_compare(slice1, slice2):
    """
    Compares two v_slices for equality
    :param slice1: A v_slice
    :param slice2: A v_slice
    :return: True if equal, False if not equal
    """
    equal = True
    if len(slice1.pitches) != len(slice2.pitches):
        equal = False
    else:
        for p in slice1.pitches:
            if p not in slice2.pitches:
                equal = False
                break
    return equal


def write_report(file, slices, general_file, general_file_command, section_name):
    """
    Writes a report to CSV
    :param file: A file name (and path if necessary)
    :param slices: A list of slices
    """

    max_density = 0  # The maximum number of pitches in a chord (may be greater than PS)
    counter = 0  # The number of values of UNS, LNS, PS, and INS to average
    ps_avg = 0  # The PS average
    uns_avg = 0  # The UNS average
    lns_avg = 0  # The LNS average
    ins_avg = 0  # The INS average
    MT_avg = 0  # The MT average
    mT_avg = 0  # The mT average
    max_ps = 0  # The PS max
    min_ps = 0  # The PS min
    max_ins = 0  # The INS max
    min_ins = 0  # The INS min
    max_uns = 0  # The UNS max
    min_uns = 0  # The UNS min
    max_lns = 0  # The LNS max
    min_lns = 0  # The LNS min
    max_MT = 0  # The MT max
    min_MT = 0  # The MT min
    max_meant = 0  # The mT max
    min_meant = 0  # The mT min

    lp = None  # The lowest pitch in the piece
    hp = None  # The highest pitch in the piece
    lps = 0  # The cardinality of LPS

    # Calculate averages, lowest pitch, highest pitch, and #(LPS)
    for s in slices:
        if s.cardinality > max_density:
            max_density = s.cardinality
        if s._ps != None:
            if s.ps > 0:
                counter += 1
                ps_avg += s.ps
                uns_avg += s.uns
                lns_avg += s.lns
                ins_avg += s.ins
                MT_avg += s.mt
                mT_avg += s.meant
        if s.lower_bound is not None and s.upper_bound is not None:
            if lp is None:
                lp = s.lower_bound
            if hp is None:
                hp = s.upper_bound
            if lp > s.lower_bound:
                lp = s.lower_bound
            if hp < s.upper_bound:
                hp = s.upper_bound
            lps = hp - lp + 1

    # Set mins to #(LPS)
    min_ps = lps
    min_uns = lps
    min_lns = lps
    min_ins = lps
    min_MT = lps
    min_meant = lps

    # Get mins and maxes
    for s in slices:
        if s.ps != None:
            if s.ps < min_ps:
                min_ps = s.ps
            if s.ps > max_ps:
                max_ps = s.ps
            if s.ps > 0:
                if s.uns < min_uns:
                    min_uns = s.uns
                if s.lns < min_lns:
                    min_lns = s.lns
                if s.ins < min_ins:
                    min_ins = s.ins
                if s.mt < min_MT:
                    min_MT = s.mt
                if s.uns > max_uns:
                    max_uns = s.uns
                if s.lns > max_lns:
                    max_lns = s.lns
                if s.ins > max_ins:
                    max_ins = s.ins
                if s.mt > max_MT:
                    max_MT = s.mt
                if s.meant is not None:
                    if s.meant < min_meant:
                        min_meant = s.meant
                    if s.meant > max_meant:
                        max_meant = s.meant

    if counter > 0:
        ps_avg /= counter
        uns_avg /= counter
        lns_avg /= counter
        ins_avg /= counter
        MT_avg /= counter
        mT_avg /= counter

    # Write the general report
    if general_file != "":
        with open(general_file, general_file_command) as general:
            if general_file_command == "w":
                # Write column headings
                general.write("Section,LPS,P_U,P_L,PS avg,PS min,PS max,UNS avg,UNS min,UNS max," + \
                              "LNS avg,LNS min,LNS max,INS avg,INS min,INS max,MT avg,MT min,MT max,mT avg,mT min,mT max\n")

            # Write general info
            general.write(section_name + ",")
            general.write(str(lps) + "," + str(hp) + "," + str(lp) + "," + str(ps_avg) + "," + \
                          str(min_ps) + "," + str(max_ps) + "," + str(uns_avg) + "," + str(min_uns) + "," + \
                          str(max_uns) + "," + str(lns_avg) + "," + str(min_lns) + "," + str(max_lns) + "," + \
                          str(ins_avg) + "," + str(min_ins) + "," + str(max_ins) + "," + str(MT_avg) + "," + \
                          str(min_MT) + "," + str(max_MT) + "," + str(mT_avg) + "," + str(min_meant) + "," + \
                          str(max_meant) + "\n")

    # Write a report to file
    with open(file, "w") as results:
        if len(slices) > 0:
            # Track the onset position in seconds
            position = 0

            # Output column headings
            line = "Measure #,Start Time (seconds),Duration (seconds),Quarter duration,Chord cardinality," + \
                   "PS,Match,NS,UNS,INS,LNS,MT,mT,Forte name,Prime form,pcset,ipseg,ipstdev"
            for i in range(max_density):
                line += ",Pitch " + str(i + 1)
            for i in range(max_ps):
                line += ",Pn_" + str(i + 1)
            line += "\n"
            results.write(line)

            # Output each slice
            for item in slices:
                c = item.chord
                s = item.get_ipseg_stdev()
                line = str(item.measure)
                line += "," + str(float(position))
                line += "," + str(float(item.duration))
                line += ",\'" + str(item.quarter_duration)
                line += "," + str(item.cardinality)
                line += "," + str(item.ps)
                if item.cardinality == item.ps:
                    line += ",TRUE"
                else:
                    line += ",FALSE"
                if item.ns is not None:
                    line += "," + str(item.ns)
                else:
                    line += ",N/A"
                if item.uns is not None:
                    line += "," + str(item.uns)
                else:
                    line += ",N/A"
                line += "," + str(item.ins)
                if item.lns is not None:
                    line += "," + str(item.lns)
                else:
                    line += ",N/A"
                if item.mt is not None:
                    line += "," + str(item.mt)
                else:
                    line += ",N/A"
                if item.meant is not None:
                    line += "," + str(item.meant)
                else:
                    line += ",N/A"
                if item.forte is not None:
                    line += ",\"\'" + str(item.forte) + "\""
                else:
                    line += ",N/A"
                if item.prime_form is not None:
                    line += ",\"" + str(item.prime_form) + "\""
                else:
                    line += ",N/A"
                if item.pcset is not None:
                    line += ",\"" + str(item.pcset) + "\""
                else:
                    line += ",N/A"
                line += "," + item.get_ipseg_string()
                if s is not None:
                    line += ",\"" + str(s) + "\""
                else:
                    line += ",N/A"
                for i in range(max_density):
                    if i < len(c.pitches):
                        line += "," + str(c.pitches[i].nameWithOctave)
                    else:
                        line += ","
                for i in range(max_ps):
                    if i < len(item.pitches_sorted):
                        line += "," + str(item.pitches_sorted[i])
                    else:
                        line += ","
                line += "\n"
                results.write(line)
                position += item.duration

        # If we have no slices to write, just write a newline
        else:
            results.write("\n")
