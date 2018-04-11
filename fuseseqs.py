import difflib
import logging as log
from seqlogging import setup_logging

setup_logging()
logger = log.getLogger(__name__)


# Considers only those sequences not included in other sequences of this same set.
def clean_patterns(S):
    Sfinal = []
    for i in range(0, len(S)):
        elements_containing_string = [index for index, item in enumerate(S) if S[i] in item]
        if len(elements_containing_string) is 1:
            Sfinal.append(S[i])
    return Sfinal


def mergeable(s1, s2):
    d = difflib.SequenceMatcher(None, s1, s2)
    match = max(d.get_matching_blocks(), key=lambda x: x[2])
    p1_start, p2_start, length = match
    if p1_start+length == len(s1) and p2_start is 0:
        log.debug("  + subsequences %s and %s are mergeable", s1, s2)
        log.debug("  + Merged as: {}".format(s1+s2[length:]))
        return [(s1+s2[length:])]
    else:
        log.debug("  - subsequences %s and %s are NOT mergeable", s1, s2)
        return None


def merge_patterns(patterns, index):
    if index >= len(patterns):
        return patterns
    m = None
    for i in range(len(patterns)):
        log.debug("Trying with index: %d", i)
        if i is not index:
            m = mergeable(patterns[index], patterns[i])
            if m is not None:
                log.debug("  + pat_list: %s", patterns)
                log.debug("  + merging.: %s (pos: %d) and %s (pos: %d)", patterns[index], index, patterns[i], i)
                first_to_be_removed = patterns[index]
                second_to_be_removed = patterns[i]
                patterns.remove(first_to_be_removed)
                patterns.remove(second_to_be_removed)
                break
            else:
                log.debug("- index %d is not mergeable!", i)
        else:
            log.debug("- ops! skipping index %d", i)
    if m is None:
        log.debug("- no more matches found.")
        return patterns
    else:
        if index+1 < len(patterns):
            log.debug("  - recursive call")
            merge_patterns(patterns+m, index+1)
        else:
            log.debug("finishing!")
            log.debug("  - patterns: %s", patterns)
            log.debug("  - m.......: %s", m)
            return patterns+m


def fuse_patterns(original_pattern_list):
    patterns_are_merging = False
    i = 0
    condition = (i < len(original_pattern_list)) or (patterns_are_merging is True)
    while condition:
        merged_patterns = merge_patterns(original_pattern_list, i)
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
        log.debug("\nNEW LIST: %s", merged_patterns)
        log.debug("----")
    log.debug("\n----\nFINAL LIST: %s", original_pattern_list)
    return original_pattern_list
