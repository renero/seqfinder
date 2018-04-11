import logging
from seqlogging import setup_logging
import collections
import ahocorasick

setup_logging()
logger = logging.getLogger(__name__)


def gen_patterns(sequence, max_wildcards, gene_length, patterns_set):
    def insert_wildcard(sequence, num_wildcards, max_wildcards, position, initial_position, cycled,
                        gene_length, patterns_set):
        # print('insert_wildcard(', sequence, ', num_wildcards=', num_wildcards,
        #       ', max_wildcards=', max_wildcards, ', position=', position,
        #       ', initial_position=', initial_position, ', cycled=', cycled,
        #       ', gene_length=', gene_length, ')', sep='')
        # print('    set: {', patterns_set, '}')

        if position >= len(sequence) and num_wildcards < max_wildcards:
            # print('  ↩︎ Cycling')
            cycled = True
            position = 0
        if position >= len(sequence) or num_wildcards >= max_wildcards or (position is initial_position and cycled):
            # print('  ✕ Out of bounds OR max num of wildcards reached OR cycled'); print()
            return

        sequence = list(sequence)

        if num_wildcards is max_wildcards:
            # print('  ⇥ Max nr. of wildcards reached', sep='')
            return sequence
        else:
            # print('  num_wildcards(', num_wildcards, ') != max_wildcards(', max_wildcards, ')', sep='')
            if sequence[position] is not '?':
                candidate = list(sequence)
                original_sequence = list(sequence)
                candidate[position:position + gene_length] = '?' * gene_length
                # print('  sequence candidate: ', ''.join(sequence), ' ▷ ', ''.join(candidate), sep='')
                sequence = candidate
                num_wildcards += 1
                if ''.join(candidate) not in patterns_set and num_wildcards is max_wildcards:
                    # print('  👍 candidate accepted: ', ''.join(sequence), sep='')
                    patterns_set.add(''.join(candidate))
                else:
                    if ''.join(candidate) in patterns_set:
                        # print('  ✕ candidate already in set: ', ''.join(candidate), sep='')
                        sequence = list(original_sequence)
                        num_wildcards -= 1
            # else:
                # print('  ⌧ Skipping')
            if position is initial_position and cycled:
                return
            else:
                return insert_wildcard(''.join(sequence), num_wildcards, max_wildcards,
                                       position + gene_length, initial_position, cycled,
                                       gene_length, patterns_set)

    patterns_set.add(sequence) # Add the initial sequence without wildcards to the set.
    for shift in range(0, len(sequence), gene_length):
        for position in range(0, len(sequence) - shift, gene_length):
            insert_wildcard(sequence, 0, max_wildcards, position + shift, position + shift, False,
                                      gene_length, patterns_set)
        # print('\n\n----\n\n')

    return patterns_set


def find_patterns(sequences, regex_length=0, gene_length=1):
    """
    Find pattern in a list of sequences. Return the winner prototype, max score and ration between
    number of matches and total number of sequences in cluster.
    pattern_length specifies the maximum number of widlcards used on each pattern. If None, the max will
    be the length of the sequence - 1.
    """
    logger.debug('Processing sequences: %s', ','.join(sequences))
    A = ahocorasick.Automaton()
    for index, word in enumerate(sequences):
        A.add_word(word, (index, word))
    max_score = -1
    max_rate = -1
    prototype = '?' * len(sequences[0])
    for seq in sequences:
        logger.debug("processing sequence: %s", seq)
        patterns_set = gen_patterns(seq, max_wildcards=regex_length, gene_length=gene_length, patterns_set=set())
        for pattern in patterns_set:
            logger.debug("  processing pattern: %s", pattern)
            score = len(list(A.keys(pattern, '?', ahocorasick.MATCH_EXACT_LENGTH)))
            rate = score / len(sequences)
            logger.debug("    score,count,rate: [%.2f,%d,%02f] .. <%.02f,%d,%.02f>",
                score, pattern.count('?'), rate, max_score, prototype.count('?'), max_rate)
            if ((score > max_score) and (rate > max_rate)) or (
                            (score == max_score) and (rate == max_rate) and (
                                pattern.count('?') < prototype.count('?'))):
                max_score = score
                prototype = pattern
                max_rate = rate
                logger.debug("      ** new max score: %.02f", max_score)
                logger.debug("      ** new prototype: %s", prototype)
                logger.debug("      ** new rate.....: %.02f", rate)
    logger.debug("WINNER ** %s", prototype)

    Cluster = collections.namedtuple('Cluster', ['prototype', 'sequences', 'score', 'rate'])
    return Cluster(prototype, sequences, max_score, max_rate)
