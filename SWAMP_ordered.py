#!/usr/bin/env python
'''
A little modification of SWAMP: write out sequences in the result philip file
in the same order as it were read.
               ____  _   _   _  ____    _______  ____
              / ___)| | | | | |/  _  \ / _   _ \|  _ \ 
              \___ \| | | | | |  (_)  | | | | | | |_) )
               ___) | |_| |_| |  ___  | | | | | |  __/
              (____/ \_______/|_|   |_|_| |_| |_|_|

         Sliding Window Alignment Masker for PAML - 31-03-14

SWAMP analyses multiple sequence alignments in a phylogenetic context,
looking for regions of higher than expected non-synonymous substitutions
along a branch, over a short sequence window.  If a user defined
threshold is exceeded then the window of sequence is masked to prevent
its inclusion in downstream evolutionary analyses.  This masking
approach removes sequence data that violates the assumptions of the
phylogenetic models implemented in the software package PAML that could
otherwise give a false signal of positive selection.

Copyright (C) 2015 Peter Harrison

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
#==============================================================================
import sys
import os
import argparse
from collections import OrderedDict
#==============================================================================
#Main==========================================================================
#==============================================================================


def main():
    '''Main runs only when the file itself is executed.
    If a part of the script is imported, main will not be executed.'''

    parser = argparse.ArgumentParser(description="SWAMP - Sliding Window \
                                     Alignment Masker for PAML - Version \
                                     31-03-14", epilog="Copyright (C) 2014 \
                                     Peter Harrison", add_help=False)

    # Required parameters
    parameters = parser.add_argument_group('Required parameters')
    parameters.add_argument("-i", "--infolder", type=str,
                            help="An input folder containing .phy and rst \
                            files. This folder can contain multiple \
                            subfolders each with .phy and rst files.")
    parameters.add_argument("-b", "--branchnames", type=str,
                            help="A file listing which branches to analyse \
                            and which sequences to mask, see documentation \
                            for file format in README.md.")
    parameters.add_argument("-t", "--threshold", type=int,
                            help="A threshold integer of the number of \
                            non-synonymous changes at and above which the \
                            window will be masked.")
    parameters.add_argument("-w", "--windowsize", type=int,
                            help="The window size (in codons) for the sliding \
                            window scan.")

    # Optional arguments
    options = parser.add_argument_group('Optional arguments')
    options.add_argument("-h", "--help", action="help",
                         help="Show help and exit.")
    options.add_argument("-m", "--minseqlength", type=int, default=33,
                         help="The required minimum number of informative \
                         codons in each sequence of the alignment \
                         post-masking, a warning will be provided if sequences\
                          are less than this length.")
    options.add_argument("-s", "--interscan", action="store_true",
                         help="Activates interscan masking as desribed in the \
                         documentation, README.md.")
    options.add_argument("-p", "--print-alignment", type=str,
                         help="Prints out a summary of the given alignment \
                         file.  No other option can be used at the same time.")

    # Check if the user supplied any arguments. If not, print help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Shorten the access to command line arguments.
    args = parser.parse_args()

    if args.print_alignment:
        print_alignment_summary(args.print_alignment)
        sys.exit(0)
    elif args.infolder is None:
        print "Parameter -i [--infolder] is required!"
        sys.exit(1)
    elif args.branchnames is None:
        print "Parameter -b [--branchnames] is required!"
        sys.exit(1)
    elif args.threshold is None:
        print "Parameter -t [--threshold] is required!"
        sys.exit(1)
    elif args.windowsize is None:
        print "Parameter -w [--windowsize] is required!"
        sys.exit(1)

    # Print user provided folder and SWAMP user options
    print "Scanning folder:\t", args.infolder
    params = (args.threshold, args.windowsize, args.minseqlength)
    print "With selected parameters:"
    print "\tThreshold:\t%s\n\tWindowsize:\t%s\n\tMinimum length:\t%s" % params
    if args.interscan:
        print "\tInterscan enabled"

    # Get all files recursively.
    file_list = list_folder_filtered(args.infolder, ".phy")

    # Ignore any already-masked files.
    file_list = [x for x in file_list if not x.endswith("_masked.phy")]

    print "Found %s phylip files to scan with SWAMP." % len(file_list)

    # Obtain sequences to be masked and branches to analyse from file
    branchcodes = read_branchcodes(args.branchnames)

    # Sliding window scan all found phylip alignments
    sliding_window_scan(file_list, args.threshold, args.windowsize,
                        args.interscan, branchcodes,
                        args.minseqlength)
#==============================================================================
#Functions=====================================================================
#==============================================================================


def print_alignment_summary(infile):
    '''Prints out a summary of the given alignment file'''
    seq_dict = read_phylip(infile)

    aln_length = len(seq_dict.values()[0])

    print "%15s  Alignment length: %d codons" % ('', aln_length)

    total_masked_sites = 0
    for species in seq_dict:
        seq = seq_dict[species]
        codons = get_codons(seq)

        n_nongap_sites = 0
        n_masked_sites = 0
        for codon in codons:
            if codon != '---':
                n_nongap_sites += 1

            if codon == 'NNN':
                n_masked_sites += 1
                total_masked_sites += 1

        pct_masked = (float(n_masked_sites) / float(n_nongap_sites)) * 100
        print ("%15s  %.50s... %5d codons %4d masked (%.0f%%)" % (species, seq,
               n_nongap_sites, n_masked_sites, pct_masked))

    print "%15s  Total masked codons: %d" % ('', total_masked_sites)


def list_folder_filtered(infolder, extension):
    '''Recursivly returns files with given extension from the specified folder
    and subfolders'''
    try:
        file_list = []

        # Walk directory tree
        for path, subdirs, files in os.walk(infolder):
            for name in files:
                if name.endswith(extension):  # Check file ends with extension
                    phyfile = os.path.join(path, name)  # Get full path
                    file_list.append(phyfile)  # Add path to file list

        # If no phylip files are found print error
        if not file_list:
            print "!----ERROR----!"
            print "%s doesn't exist or does not contain .phy files!" % infolder
            sys.exit(1)

        return file_list

    except KeyboardInterrupt:
        sys.exit(1)


def read_branchcodes(infile):
    '''Parses branch information file for sequences to mask and branches to use
    for non-synonymous information, see documentation for file format.'''
    try:
        branchcodes = {}

        with open(infile) as readfile:
            for line in readfile:
                parts = line.split()

                # Check line has branch and sequence names
                if len(parts) == 2:
                    seqparts = parts[1].split(",")  # Obtain sequences to mask
                    branchcodes[parts[0]] = tuple(seqparts)  # Store branches

                # Raise error for lines that don't have two parts
                else:
                    raise ValueError(("The following line in %s is not "
                                      "formatted correctly. Check the "
                                      "documentation:\n%s" % (infile, line)))
        return branchcodes

    except IOError:
        print "File %s does not exist!" % infile
        raise

    except KeyboardInterrupt:
        sys.exit(1)


def sliding_window_scan(file_list, threshold, windowsize,
                        interscan, branchcodes, minseqlength=None):
    '''Sliding window scan and mask a list of files'''

    # Counter for total masked and short seqs
    masked_file_count, files_with_short_seqs_count = 0, 0

    for infile in file_list:
        result = sliding_window_scan_file(infile, threshold, windowsize,
                                          interscan, branchcodes, minseqlength)

        if result['masked_codon_count'] > 0:
            masked_file_count += 1

        if result['short_seq_count'] > 0:
            files_with_short_seqs_count += 1
            if __name__ == "__main__":
                sinfo = (result['short_seq_count'], minseqlength)
                print "    !! %s sequence(s) contained fewer than %s" % sinfo
                print "    informative codons after masking!"

    if __name__ == "__main__":
        filesinfo = (masked_file_count, len(file_list))
        print "SWAMP has masked %s out of %s PAML phylip files" % filesinfo
        if minseqlength:
            pin = (files_with_short_seqs_count, minseqlength)
            print "%s files have sequences less than %s codons in length" % pin


def sliding_window_scan_file(infile, threshold, windowsize, interscan,
                             branchcodes, minseqlength=None):
    '''Sliding window scan and mask a single file. Returns a 'result'
    dict with info resulting from the scan & masking process.'''

    result = {
        'masked_codon_count': 0,
        'masked_column_count': 0,
        'short_seq_count': 0,
        'masked_dict': None,
        'seq_dict': None
    }

    seqdict = read_phylip(infile)  # Create sequences dictionary
    branch_error_check(branchcodes, seqdict)  # Check validity of branchcodes

    # Extract branch information from rst file
    branches = read_rst(infile)

    # Scan codons for those that need masking
    codons_tomask = (scan(branches, seqdict, threshold, windowsize,
                          branchcodes))

    # Check if any codons were needing to be masked
    if len(codons_tomask) > 0:

        # If interscan option was selected
        if interscan:
            codons_tomask = run_interscan(
                codons_tomask, seqdict, branchcodes)

    # Mask the codons to get a new seq_dict
    masked_dict = mask_codons(seqdict, codons_tomask)
    # Count number of codons with >0 species masked out.
    result['masked_column_count'] = len(codons_tomask.keys())
    result['masked_codon_count'] = count_masked_codons(masked_dict)

    # Store the # of sequences with less than the desired minimum number of
    # codons
    if minseqlength is not None:
        short_seq_count = count_short_seqs(masked_dict, minseqlength)
        result['short_seq_count'] = short_seq_count

    if __name__ == "__main__":
        n_species = len(seqdict.keys())
        n_columns = len(seqdict.values()[0]) / 3
        spinfo = (infile, n_species, n_columns)
        print "  > %s  (%d species, %d columns)" % spinfo
        minfo = (result['masked_codon_count'], result['masked_column_count'])
        print "    masked %d codons from %d columns" % minfo

    # Print out masked alignment.
    print_masked_phyfile(infile, masked_dict)

    result['masked_dict'] = masked_dict
    result['seq_dict'] = seqdict

    return result


def count_masked_codons(seq_dict):
    '''Count number of masked codons'''
    masked_codon_count = 0
    for species in seq_dict:
        seq = seq_dict[species]
        codons = get_codons(seq)
        for codon in codons:
            if codon == 'NNN':
                masked_codon_count += 1

    return masked_codon_count


def read_phylip(infile):
    '''Read in phylip sequences'''
    try:
        with open(infile) as readfile:
            # sequences = {}
            sequences = OrderedDict()
            # Get header
            header = next(readfile).split()

            # Check the supplied PHYLIP file is valid
            if len(header) != 2:
                raise ValueError("Cannot detect header in %s" % infile)
            else:
                expected_number = int(header[0])
                seq_len = int(header[1])

            # Parse the phylip by tracking a 'reading_seq' state
            # variable, which flips to False when we expect to see
            # another sequence ID.
            cur_name = None
            reading_seq = False

            for line in readfile:
                line = line.rstrip()

                if not reading_seq:
                    # Get here if we're expecting a sequence ID. Take
                    # the current line as the ID and create an entry
                    # in the sequences dict with an empty string.
                    cur_name = line
                    sequences[cur_name] = ""
                    reading_seq = True
                else:
                    # Get here if we're expecting sequence data.
                    sequences[cur_name] += line.upper()
                    cur_seq_len = len(sequences[cur_name])
                    if cur_seq_len == seq_len:
                        # Once the sequence has reached the
                        # header-defined length, we know to expect a
                        # sequence ID on the next line.
                        reading_seq = False

            # Check if sequences are correct
            if len(sequences) != expected_number:
                err = (expected_number, len(sequences), infile)
                raise ValueError(
                    "Expected %s sequences and have loaded %s in %s" % err)

            # Check sequence length and frame
            for seq in sequences:

                if len(sequences[seq]) > seq_len:
                    raise ValueError(
                        "Sequence lengths are not all equal for: %s" % infile)

                if len(sequences[seq]) == 0:
                    raise ValueError("Sequence must be longer than 0!")

                if len(sequences[seq]) % 3 != 0:
                    raise ValueError("%s sequence not divisible by 3" % seq)
        return sequences

    except IOError:
        print "File %s does not exist!" % infile
        raise

    except KeyboardInterrupt:
        sys.exit(1)


def branch_error_check(branchcodes, seqdict):
    '''Check for errors in branchcodes file and loaded branches'''
    # If there are no valid branchcodes
    if not branchcodes:
        raise ValueError("No branchnames found in --branchnames file")

    # Check for errors in supplied branchcodes
    for branch in branchcodes:
        if ".." not in branch:
            raise ValueError("Check branchnames file. Branch %s is invalid" %
                             branch)

        # Given list of species from branch code...
        species_list = branchcodes[branch]
        for species in species_list:
            # Raise an error if any species isn't in the alignment dict.
            if species not in seqdict:
                raise ValueError("Species %s not in alignment file!" % species)


def read_rst(infile):
    '''Parse rst file for non-synonymous substitutions'''

    # Get path for rst file, should be in same folder as phylip file
    rstinfile = "/".join(infile.split("/")[:-1])+"/rst"

    try:
        with open(rstinfile) as readfile:
            branches = {}
            readin = False  # Save or not save line toggle
            no_branch_found = True  # No branch file check

            for line in readfile:
                line = line.rstrip()

                # Find branch line
                if line.startswith("Branch"):
                    parts = line.split()
                    branch = parts[2]
                    branches[branch] = []
                    readin = True
                    no_branch_found = False

                # End of branch section
                elif line.startswith("List"):
                    readin = False

                elif readin is True:
                    parts = line.split()
                    if len(parts) > 1:  # If line not empty
                        if parts[2] != parts[6]:  # If it is non-synonymous
                            codon = int(parts[0])  # Get non-synonymous codon
                            branches[branch].append(codon)

            if no_branch_found:
                print "Check that PAML ran correctly!"
                raise ValueError(
                    "File %s is missing branch information!" % rstinfile)

        return branches
    except IOError:
        print "File %s does not exist!" % rstinfile
        raise
    except KeyboardInterrupt:
        sys.exit(1)


def scan(branches, seqdict, threshold, windowsize, branchcodes):
    '''Conducts sliding window scan on user supplied sequences and branches.'''
    # Dictionary to store codons that need masking
    codons_tomask = {}

    # For each of the branches to be analysed
    for branch in branchcodes:
        if branch not in branches:
            print "!----ERROR----!"
            print "Check branchnames file as branch not valid: ", branch
            sys.exit(1)
        codons = branches[branch]
        windowstart = 0
        windowend = windowsize - 1

        # For each window
        while windowend <= len(seqdict.values()[0]) / 3:
            count = 0
            window = range(windowstart, windowend)  # Get codons in window

            for codon in codons:
                if codon in window:  # Count substitions within window
                    count += 1

            # If threshold equaled or exceeded
            if count >= threshold:
                for codon in window:
                    if codon in codons_tomask:  # Mark window for masking
                        for seq in branchcodes[branch]:
                            if seq not in codons_tomask[codon]:
                                codons_tomask[codon].append(seq)
                    else:
                        codons_tomask[codon] = list(branchcodes[branch])

            windowstart += 1
            windowend += 1

    return codons_tomask


def run_interscan(codons_tomask, seqdict, branchcodes):
    '''Extend masking around masked sections according to length'''
    seqs_tomask = []

    # Get list of sequences that need masking
    for branch in branchcodes.values():
        for seq in branch:
            seqs_tomask.append(seq)

    seqs_tomask = list(set(seqs_tomask))
    seqlength = len(seqdict.values()[0])/3

    # Scan each sequence seperately, but only those in branchcodes file
    for seq in seqs_tomask:
        startat = 0  # Start of sequence
        scan_complete = False  # Boolean for completion

        # Only use codons specific to the sequence in question
        bs_codons_tomask = {}
        for codon in codons_tomask:
            if seq in codons_tomask[codon]:
                bs_codons_tomask[codon] = seq

        #Run until no more sections to mask
        while not scan_complete:
            bs_codons_tomask, scan_complete, startat = (ngapscan(
                bs_codons_tomask, seqlength, scan_complete, startat))

        # Scan in reverse for end section until no more changes are found
        rev_scan_complete = False
        while not rev_scan_complete:  # Scan in reverse for end section
            bs_codons_tomask, changemade = (ngapscan_reverse(
                bs_codons_tomask, seqlength))

            # If change made scan in forward direction again
            if changemade:
                startat = 0  # Start from begining again
                scan_complete = False

                while not scan_complete:
                    bs_codons_tomask, scan_complete, startat = (ngapscan(
                        bs_codons_tomask, seqlength, scan_complete, startat))

            # If no more sections can be masked
            else:
                rev_scan_complete = True

        # Add new branch specific intermasking to main codons_tomask dictionary
        for codon in bs_codons_tomask:
            if codon in codons_tomask:
                if seq not in codons_tomask[codon]:
                    codons_tomask[codon].append(seq)
            else:
                codons_tomask[codon] = [seq]

    return codons_tomask


def ngapscan(codons_tomask, seqlength, scan_complete, startat):
    '''Conduct intermask scan and masking'''
    # Set starting window position
    current = startat
    start, masked1, inter, masked2 = 0, 0, 0, 0  # Set codon counters to 0
    to_mask_start, to_mask_inter = [], []  # Lists to hold codons

    # Scan for length of unmasked region at start of sequence if present
    while current <= seqlength:
        if current in codons_tomask:  # Found first masked section
            break
        else:
            start += 1
            to_mask_start.append(current)
        current += 1

    # Scan length of masked region
    while current <= seqlength:
        if current not in codons_tomask:
            break
        else:
            masked1 += 1
        current += 1

    # Scan for length of intermasked (unmasked) region
    while current <= seqlength:
        if current in codons_tomask:
            newstart = current
            break
        else:
            inter += 1
            to_mask_inter.append(current)
        current += 1

    # Scan for length of second masked region
    while current <= seqlength:
        if current not in codons_tomask:
            break
        else:
            masked2 += 1
        current += 1

    # If unmasked region at start of sequence less than twice length of
    # masked region then mask it.
    if start > 0 and start < (masked1*2):
        for codon in to_mask_start:
            codons_tomask[codon] = None
            startat = 0

    # If second masked region exists
    elif masked2 > 0:

        # If intermasked region less length than masked regions either
        # side of it then mask it
        if inter < (masked1+masked2):
            for codon in to_mask_inter:
                codons_tomask[codon] = None
                startat = 0  # Scan from start again
        else:  # Look for more intermasked regions in the same sequence
            startat = newstart
            scan_complete = True

    # No more intermasked regions to scan
    else:
        scan_complete = True

    return codons_tomask, scan_complete, startat


def ngapscan_reverse(codons_tomask, seqlength):
    '''Conduct intermask scan and masking in reverse'''
    # Start at end of sequence
    current = seqlength
    end, masked1 = 0, 0  # Set codon counters to 0
    to_mask_end = []

    # Scan for length of unmasked region at end of sequence if present
    while current >= 0:
        if current in codons_tomask:  # Found first masked section
            break
        else:
            end += 1
            to_mask_end.append(current)
        current -= 1

    # Get length of masked region
    while current <= seqlength:
        if current not in codons_tomask:
            break
        else:
            masked1 += 1
        current -= 1

    # If unmasked region at end of sequence less than twice length of masked
    #region next to it then mask it.
    if end > 0 and end < (masked1*2):
        for codon in to_mask_end:
            codons_tomask[codon] = None
        changemade = True
    else:
        changemade = False

    return codons_tomask, changemade


def print_masked_phyfile(phyinfile, masked_dict):
    '''Masks codons and prints the Phylip output file.'''

    # Create the output file and write Phylip data.
    phyoutfile = phyinfile.replace(".phy", "_masked.phy")
    print_phyfile(phyoutfile, masked_dict)

    if __name__ == "__main__":
        print "    wrote file %s" % phyoutfile


def mask_codons(seqdict, codons_tomask):
    '''Masks a seq dict given the codons_tomask datastructure, which
    is keyed by codon index with each value being a tuple of the seq
    IDs needing masking.'''
    masked_dict = OrderedDict()
    for seq in seqdict.keys():
        codons = get_codons(seqdict[seq])  # Split sequence into codons
        i = 1
        newsequence = ""

        # For each codon in sequence
        for codon in codons:
            # If codon and sequence is to be masked
            if i in codons_tomask and seq in codons_tomask[i]:
                newsequence += "NNN"  # Replace codon with N's
            else:  # If codon is not to be masked
                newsequence += codon  # Include original codon unmasked
            i += 1

        # Store sequences for length check
        masked_dict[seq] = newsequence

    return masked_dict


def print_phyfile(phyoutfile, seqdict):
    '''Prints a sequence dict in phylip format.'''
    out = create_outfile(phyoutfile)  # Create output masked phylip file

    # Print phylip header, needs number of sequences and sequence length
    print >>out, "      "+str(len(seqdict))+"   "+str(len(seqdict.values()[0]))

    # for seq_id in sorted(seqdict.keys()):
    for seq_id in seqdict.keys():
        print >>out, str(seq_id)  # Print sequence header
        seq_str = seqdict[seq_id]
        lines = format_as_phylip(seq_str)

        # Print each line of masked sequence
        for line in lines:
            print >>out, line


def count_short_seqs(masked_dict, minseqlength):
    '''Count the number of sequences with fewer than minseqlength
    codons post-masking.'''

    short_seq_count = 0

    # For each sequence
    for seq in sorted(masked_dict.keys()):
        masked_codons = {}
        codons = get_codons(masked_dict[seq])  # Get codons
        i = 0

        while i < len(codons):
            if "N" in codons[i]:  # If codon has N charachter
                masked_codons[i] = 1
            i += 1

        # Obtain length of informative sites
        unmasked_length = len(codons)-len(masked_codons)

        if unmasked_length < minseqlength:
            short_seq_count += 1

    return short_seq_count


def create_outfile(outfile):
    '''Create outfile'''
    try:
        out = open(outfile, 'w')
        return out

    except IOError:
        print "File %s cannot be created!" % outfile
        raise

    except KeyboardInterrupt:
        sys.exit(1)


def get_codons(sequence):
    '''Split DNA sequence into codons'''
    codon_list = []

    # Divide sequence into triplets
    for i in range(0, len(sequence), 3):
        codon_list.append(sequence[i: i + 3])  # Add codons to list

    return codon_list


def format_as_phylip(line):
    '''Split sequence into 60 character lines for printing in phylip format'''
    phylip_lines = []

    for i in range(0, len(line), 60):
        phylip_lines.append(line[i: i + 60])  # Add 60 character chunks to list

    return phylip_lines

if __name__ == "__main__":
    main()
