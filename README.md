# VerticalAnalyzer
Written by Jeff Martin
Copyright (c) 2022 by Jeff Martin. All rights reserved. This code is made available under the GNU GPL v3.0. See the LICENSE file for more information.

## Description
This program uses the music21 library to parse musicXML files. It then segments the score by chord (simultaneity) and calculates many interesting statistics.

## List of data extracted from a score for each sonority
- Measure number
- Start time (in seconds) relative to the start of the piece
- Duration (in seconds)
- Quarter duration
- Cardinality of chord (includes duplicate pitches)
- Cardinality of pitch set (excludes duplicate pitches)
- Duplicates present (TRUE/FALSE)
- Set-class name (Forte, Rahn, and Carter)
- Core harmony (used in Carter's music, includes all-interval tetrachords and the all-trichord hexachord)
- Derived core harmony (after John Link)
- Derived core harmony associations
- Pcset (pitch-class set)
- Pset (pitch set)
- PSC (pitch set-class)
- PSCS (pset spacing contour segment)
- PSI (pset spacing index)
- List of pitch names and integers in the chord, from low to high

## List of data extracted from a score for each movement, as well as the entire work
- Starting time (in seconds) relative to the start of the piece
- Duration (in seconds)
- Average pset cardinality
- Average PSI
- Interval between highest and lowest pitches
- Highest pitch
- Lowest pitch
- Frequency of occurrence of each pitch-class
- Cumulative duration of each pitch-class
- Frequency of occurrence of each pitch
- Cumulative duration of each pitch
