import csv
import difflib
import logging as log
import sys


def read_results(results_file, threshold=0.0):
    """ Reads the results file (seq, #a,#b,diff,mag,f5,threshold,tg) """
    patterns = []
    with open(results_file, 'r') as input_file:
        patterns_reader = csv.reader(input_file)
        next(input_file)
        for line in patterns_reader:
            #if float(line[6]) >= threshold:
            patterns.append(line[0])
    print("Read {:d} patterns.".format(len(patterns)))
    return patterns


# Considers only those sequences not included in other sequences of this same set.
def clean_patterns(S):
    Sfinal = []
    for i in range(0, len(S)):
        elements_containing_string = [index for index, item in enumerate(S) if S[i] in item]
        if len(elements_containing_string) is 1:
            Sfinal.append(S[i])
    log.info('Cleaned list: %s', Sfinal)
    return Sfinal


def mergeable(s1, s2):
    d = difflib.SequenceMatcher(None, s1, s2)
    match = max(d.get_matching_blocks(), key=lambda x: x[2])
    p1_start, p2_start, length = match
    if p1_start+length == len(s1) and p2_start is 0:
        log.info("  + subsequences %s and %s are mergeable", s1, s2)
        log.info("  + Merged as: {}".format(s1+s2[length:]))
        return [(s1+s2[length:])]
    else:
        log.info("  - subsequences %s and %s are NOT mergeable", s1, s2)
        return None


def merge_patterns(patterns, index):
    if index >= len(patterns):
        return patterns
    m = None
    for i in range(len(patterns)):
        log.info("Trying with index: %d", i)
        if i is not index:
            m = mergeable(patterns[index], patterns[i])
            if m is not None:
                log.info("  + pat_list: %s", patterns)
                log.info("  + merging.: %s (pos: %d) and %s (pos: %d)", patterns[index], index, patterns[i], i)
                first_to_be_removed = patterns[index]
                patterns.remove(first_to_be_removed)
                if i is not index and i < len(patterns):
                    second_to_be_removed = patterns[i]
                    patterns.remove(second_to_be_removed)
                break
            else:
                log.info("- index %d is not mergeable!", i)
        else:
            log.info("- ops! skipping index %d", i)
    if m is None:
        log.info("- no more matches found.")
        return patterns
    else:
        if index+1 < len(patterns):
            log.debug("  - recursive call")
            merge_patterns(patterns+m, index+1)
        else:
            log.info("finishing!")
            log.info("  - patterns: %s", patterns)
            log.info("  - m.......: %s", m)
            return patterns+m


def fuse_patterns(original_pattern_list):
    patterns_are_merging = False
    i = 0
    condition = (i < len(original_pattern_list)) or (patterns_are_merging is True)
    while condition:
        merged_patterns = merge_patterns(original_pattern_list, i)
        log.info("ORIGINAL LIST: %s", merged_patterns)
        if merged_patterns is original_pattern_list:
            patterns_are_merging = False
        else:
            patterns_are_merging = True
            if merged_patterns is not None:
                original_pattern_list = merged_patterns
            i = 0
        i += 1
        if len(original_pattern_list) > i or patterns_are_merging is True:
            condition = True
        else:
            condition = False
        log.info("\nNEW LIST: %s", merged_patterns)
        log.info("----")
    log.info("\n----\nFINAL LIST: %s", original_pattern_list)
    return original_pattern_list


#
# main.
#
def fuse_sequences(infile, outfile):
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.WARN)
    # if outfile is '':
    #     outputfile = sys.stdout
    # else:
    #     outputfile = open(outfile, 'w')
    # writer = csv.writer(outputfile, delimiter=',')
    raw_sequences = read_results(infile)
    log.info('Raw blocks..: %s', raw_sequences)
    patterns = clean_patterns(raw_sequences)
    merged_patterns = fuse_patterns(patterns)
    cleaned_merged_patterns = set(clean_patterns(merged_patterns))
    log.info('Fused blocks: %s', cleaned_merged_patterns)
    print("\n".join(merged_patterns))


if (2 <= len(sys.argv) <= 3) is False:
    print(len(sys.argv))
    print("Usage: {} patterns_file [output_file]".format(sys.argv[0]))
    exit(1)
if len(sys.argv) is 2:
    outfile = ''
else:
    outfile = sys.argv[2]
fuse_sequences(sys.argv[1], outfile)
