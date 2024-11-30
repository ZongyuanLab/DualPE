import os
import sys
import re
import time
import logging


def main(target_name,target_sequence,pe_format,editing_window_left, editing_window_right,pbs_length_list,homology_overlap,filter_tms,Exclude_first_C,pbs_Tm_Recommended):


    ##### Initialize arguments

    # Default PBS and RTT lengths to design
    # if pbs_length_list == 0:
    #     pbs_length_list = list(range(9, 16))
    #
    #
    # # Output directory date and time stamped
    # homology_overlap = 30

    out_dir = '%s_PrimeDesign' % str(time.strftime("%y%m%d_%H.%M.%S", time.localtime()))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Initialize logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    fh = logging.FileHandler(out_dir + '/PrimeDesign.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

 ##### IUPAC code map
    iupac2bases_dict = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G', 'a': 'a', 't': 't', 'c': 'c', 'g': 'g',
                        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
                        'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACTG]',
                        'r': '[ag]', 'y': '[ct]', 's': '[gc]', 'w': '[at]', 'k': '[gt]', 'm': '[ac]', 'b': '[cgt]',
                        'd': '[agt]', 'h': '[act]', 'v': '[acg]', 'n': '[actg]',
                        '(': '(', ')': ')', '+': '+', '-': '-', '/': '/', '.': '.'}

    def iupac2bases(iupac):

        try:
            bases = iupac2bases_dict[iupac]
        except:
            logger.error('Symbol %s is not within the IUPAC nucleotide code ...' % str(iupac))
            sys.exit(1)

        return (bases)

    # GC content
    def gc_content(sequence):
        sequence = sequence.upper()
        GC_count = sequence.count('G') + sequence.count('C')
        GC_content = float(GC_count) / float(len(sequence))

        return ("%.2f" % GC_content)

    def PBS_Tm(sequence):
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')

        score = gc_count * 4 + at_count * 2
        return score

    pbs_length_list_len = len(pbs_length_list)
    def recommend_Tm(pbs_Tm, previous_Tms):
        if pbs_Tm == pbs_Tm_Recommended:
            return True
        elif pbs_Tm == pbs_Tm_Recommended + 2 and (len(previous_Tms) < pbs_length_list_len or previous_Tms[-pbs_length_list_len] != pbs_Tm_Recommended):
            return True
        else:
            return False


    # Reverse complement function
    def reverse_complement(sequence):
        sequence = sequence
        new_sequence = ''
        for base in sequence:
            if base == 'A':
                new_sequence += 'T'
            elif base == 'T':
                new_sequence += 'A'
            elif base == 'C':
                new_sequence += 'G'
            elif base == 'G':
                new_sequence += 'C'
            elif base == 'a':
                new_sequence += 't'
            elif base == 't':
                new_sequence += 'a'
            elif base == 'c':
                new_sequence += 'g'
            elif base == 'g':
                new_sequence += 'c'
            elif base == '[':
                new_sequence += ']'
            elif base == ']':
                new_sequence += '['
            elif base == '+':
                new_sequence += '+'
            elif base == '-':
                new_sequence += '-'
            elif base == '/':
                new_sequence += '/'
            elif base == '(':
                new_sequence += ')'
            elif base == ')':
                new_sequence += '('
            elif base == '.':
                new_sequence += '.'
        return (new_sequence[::-1])

# Amino acid code
    codon_dict = {
        'GGG': ['Gly', 'G', 0.25], 'GGA': ['Gly', 'G', 0.25], 'GGT': ['Gly', 'G', 0.16], 'GGC': ['Gly', 'G', 0.34],
        'GAG': ['Glu', 'E', 0.58], 'GAA': ['Glu', 'E', 0.42], 'GAT': ['Asp', 'D', 0.46], 'GAC': ['Asp', 'D', 0.54],
        'GTG': ['Val', 'V', 0.47], 'GTA': ['Val', 'V', 0.11], 'GTT': ['Val', 'V', 0.18], 'GTC': ['Val', 'V', 0.24],
        'GCG': ['Ala', 'A', 0.11], 'GCA': ['Ala', 'A', 0.23], 'GCT': ['Ala', 'A', 0.26], 'GCC': ['Ala', 'A', 0.4],
        'AGG': ['Arg', 'R', 0.2], 'AGA': ['Arg', 'R', 0.2], 'AGT': ['Ser', 'S', 0.15], 'AGC': ['Ser', 'S', 0.24],
        'AAG': ['Lys', 'K', 0.58], 'AAA': ['Lys', 'K', 0.42], 'AAT': ['Asn', 'N', 0.46], 'AAC': ['Asn', 'N', 0.54],
        'ATG': ['Met', 'M', 1], 'ATA': ['Ile', 'I', 0.16], 'ATT': ['Ile', 'I', 0.36], 'ATC': ['Ile', 'I', 0.48],
        'ACG': ['Thr', 'T', 0.12], 'ACA': ['Thr', 'T', 0.28], 'ACT': ['Thr', 'T', 0.24], 'ACC': ['Thr', 'T', 0.36],
        'TGG': ['Trp', 'W', 1], 'TGA': ['End', 'X', 0.52], 'TGT': ['Cys', 'C', 0.45], 'TGC': ['Cys', 'C', 0.55],
        'TAG': ['End', 'X', 0.2], 'TAA': ['End', 'X', 0.28], 'TAT': ['Tyr', 'Y', 0.43], 'TAC': ['Tyr', 'Y', 0.57],
        'TTG': ['Leu', 'L', 0.13], 'TTA': ['Leu', 'L', 0.07], 'TTT': ['Phe', 'F', 0.45], 'TTC': ['Phe', 'F', 0.55],
        'TCG': ['Ser', 'S', 0.06], 'TCA': ['Ser', 'S', 0.15], 'TCT': ['Ser', 'S', 0.18], 'TCC': ['Ser', 'S', 0.22],
        'CGG': ['Arg', 'R', 0.21], 'CGA': ['Arg', 'R', 0.11], 'CGT': ['Arg', 'R', 0.08], 'CGC': ['Arg', 'R', 0.19],
        'CAG': ['Gln', 'Q', 0.75], 'CAA': ['Gln', 'Q', 0.25], 'CAT': ['His', 'H', 0.41], 'CAC': ['His', 'H', 0.59],
        'CTG': ['Leu', 'L', 0.41], 'CTA': ['Leu', 'L', 0.07], 'CTT': ['Leu', 'L', 0.13], 'CTC': ['Leu', 'L', 0.2],
        'CCG': ['Pro', 'P', 0.11], 'CCA': ['Pro', 'P', 0.27], 'CCT': ['Pro', 'P', 0.28], 'CCC': ['Pro', 'P', 0.33],
    }

    # Create codon swap dictionaries
    aa2codon = {}
    for codon in codon_dict:
        if codon_dict[codon][1] not in aa2codon:
            aa2codon[codon_dict[codon][1]] = []

        aa2codon[codon_dict[codon][1]].append([codon, codon_dict[codon][2]])

    for codon in aa2codon:
        aa2codon[codon] = sorted(aa2codon[codon], key=lambda x: x[1], reverse=True)

    # sorted aa2codon {'K': [['AAG', 0.58], ['AAA', 0.42]],'N': [['AAC', 0.54], ['AAT', 0.46]]}

    codon_swap_0 = {}
    codon_swap_1_1 = {}
    codon_swap_1_2 = {}
    codon_swap_2 = {}
    for codon in codon_dict:

        codon_swap_0[codon] = []
        codon_swap_1_1[codon] = []
        codon_swap_1_2[codon] = []
        codon_swap_2[codon] = []

        for other_codon in aa2codon[codon_dict[codon][1]]:

            # Check if PAM disrupted with silent mutations
            if codon[1:] != other_codon[0][1:]:
                codon_swap_0[codon].append(other_codon)

            if codon[2:] != other_codon[0][2:]:
                codon_swap_1_1[codon].append(other_codon)

            if codon[:1] != other_codon[0][:1]:
                codon_swap_1_2[codon].append(other_codon)

            if codon[:2] != other_codon[0][:2]:
                codon_swap_2[codon].append(other_codon)

    for codon in codon_dict:
        codon_swap_0[codon] = sorted(codon_swap_0[codon], key=lambda x: x[1], reverse=True)
        codon_swap_1_1[codon] = sorted(codon_swap_1_1[codon], key=lambda x: x[1], reverse=True)
        codon_swap_1_2[codon] = sorted(codon_swap_1_2[codon], key=lambda x: x[1], reverse=True)
        codon_swap_2[codon] = sorted(codon_swap_2[codon], key=lambda x: x[1], reverse=True)

##### Extract reference and edited sequence information
    def process_sequence(input_sequence):
        # Check formatting is correct
        format_check = ''
        for i in input_sequence:
            if i == '(':
                format_check += '('
            elif i == ')':
                format_check += ')'
            elif i == '/':
                format_check += '/'
            elif i == '+':
                format_check += '+'
            elif i == '-':
                format_check += '-'
            elif i == '.':
                format_check += '.'

        # Check composition of input sequence
        if len(input_sequence) != sum(
                [1 if x in ['A', 'T', 'C', 'G', '(', ')', '+', '-', '/', '.'] else 0 for x in input_sequence.upper()]):
            logger.error(
                'Input sequence %s contains a character not in the following list: A,T,C,G,(,),+,-,/ ...' % str(
                    input_sequence))
            sys.exit(1)

        # Check formatting
        if format_check.count('(') == format_check.count(')') and format_check.count(
                '(') > 0:  # Left and right parantheses equal
            if '((' not in format_check:  # Checks both directions for nested parantheses
                if '()' not in format_check:  # Checks for empty annotations
                    if sum([1 if x in format_check else 0 for x in
                            ['++', '--', '//', '+-', '+/', '-+', '-/', '/+', '/-', '/(', '+(', '-(', ')/', ')+',
                             ')-']]) == 0:
                        pass
                    else:
                        logger.error(
                            'Input sequence %s has more than one edit annotation per parantheses set (i.e. //,  +- , -/, etc.) ...' % str(
                                input_sequence))
                        sys.exit(1)
                else:
                    logger.error(
                        'Input sequence %s has empty parantheses without an edit annotation (i.e. /,  + , -) ...' % str(
                            input_sequence))
                    sys.exit(1)
            else:
                logger.error('Input sequence %s has nested parantheses which is not allowed ...' % str(input_sequence))
                sys.exit(1)
        else:
            logger.error('Input sequence %s does not have full sets of parantheses ...' % str(input_sequence))
            sys.exit(1)

        # Create mapping between input format and reference and edit sequence
        editformat2sequence = {}
        edits = re.findall('\(.*?\)', input_sequence)
        for edit in edits:
            if '/' in edit:  # editformat2sequence = {'(A/T)': ['A', 'T']}
                editformat2sequence[edit] = [edit.split('/')[0].replace('(', ''), edit.split('/')[1].replace(')', '')]
            elif '+' in edit:  # editformat2sequence = {'(A+T)': ['', 'T']}
                editformat2sequence[edit] = ['', edit.split('+')[1].replace(')', '')]
            elif '-' in edit:  # editformat2sequence = {'(A-T)': ['T', '']}
                editformat2sequence[edit] = [edit.split('-')[1].replace(')', ''), '']

        # Create mapping between edit number and reference and edit sequence
        editformat2sequence = {}
        editnumber2sequence = {}
        edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
        edit_counter = 1
        for edit_idx in edit_idxs:
            edit = input_sequence[edit_idx[0]:edit_idx[1]]

            if '/' in edit:
                editformat2sequence[edit] = [edit.split('/')[0].replace('(', ''),
                                             edit.split('/')[1].replace(')', '').lower(), edit_counter]
                editnumber2sequence[edit_counter] = [edit.split('/')[0].replace('(', ''),
                                                     edit.split('/')[1].replace(')', '').lower()]

            elif '+' in edit:
                editformat2sequence[edit] = ['', edit.split('+')[1].replace(')', '').lower(), edit_counter]
                editnumber2sequence[edit_counter] = ['', edit.split('+')[1].replace(')', '').lower()]

            elif '-' in edit:
                editformat2sequence[edit] = [edit.split('-')[1].replace(')', ''), '', edit_counter]
                editnumber2sequence[edit_counter] = [edit.split('-')[1].replace(')', ''), '']

            edit_counter += 1

        edit_start = min([i.start() for i in re.finditer('\(', input_sequence)])
        edit_stop = max([i.start() for i in re.finditer('\)', input_sequence)])

        edit_span_sequence_w_ref = input_sequence[edit_start:edit_stop + 1]
        edit_span_sequence_w_edit = input_sequence[edit_start:edit_stop + 1]


        for edit in editformat2sequence:
            edit_span_sequence_w_ref = edit_span_sequence_w_ref.replace(edit, editformat2sequence[edit][0])
            edit_span_sequence_w_edit = edit_span_sequence_w_edit.replace(edit, editformat2sequence[edit][1])

        edit_start_in_ref = re.search('\(', input_sequence).start()
        edit_stop_in_ref_rev = re.search('\)', input_sequence[::-1]).start()

        edit_span_length_w_ref = len(edit_span_sequence_w_ref)

        edit_span_length_w_edit = len(edit_span_sequence_w_edit)
        is_shorter_than_30 = len(edit_span_sequence_w_edit) < 30
        reference_sequence = input_sequence
        edit_sequence = input_sequence
        editnumber_sequence = input_sequence
        for edit in editformat2sequence:
            reference_sequence = reference_sequence.replace(edit, editformat2sequence[edit][0])
            edit_sequence = edit_sequence.replace(edit, editformat2sequence[edit][1])
            editnumber_sequence = editnumber_sequence.replace(edit, str(editformat2sequence[edit][2]))

        return (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence,
                edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev, is_shorter_than_30)

    ##### Dictionary for to organize different DNA targets
    target_design = {}



    target_sequence = target_sequence.upper()
    editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev, is_shorter_than_30 = process_sequence(
        target_sequence)

    # Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
    target_design[target_name] = {'target_sequence': target_sequence,
                                  'editformat2sequence': editformat2sequence,
                                  'editnumber2sequence': editnumber2sequence,
                                  'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence,
                                  'editnumber_sequence': editnumber_sequence,
                                  'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit],
                                  'edit_start_in_ref': edit_start_in_ref,
                                  'edit_stop_in_ref_rev': edit_stop_in_ref_rev,
                                  'is_shorter_than_30': is_shorter_than_30,
                                  'pegRNA': {'+': [], '-': []},
                                  }


    if len(target_design) == 0:
        logger.error(
            'Input file %s does not have any entries. Make sure a column header is included (target_name,target_sequence) ...'
                )
        sys.exit(1)

    ##### Find cut index and reformat PE format parameter
    # 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
    if (pe_format.count('[') + pe_format.count(']')) == 2:

        if pe_format.count('/') == 1:

            # Find indices but shift when removing annotations
            cut_idx = re.search('/', pe_format).start()
            pam_start_idx = re.search('\[', pe_format).start()
            pam_end_idx = re.search('\]', pe_format).start()

            # Find pam and total PE format search length
            pam_length = pam_end_idx - pam_start_idx - 1
            pe_format_length = len(pe_format) - 3

            # Check if cut site is left of PAM
            if cut_idx < pam_start_idx:

                # Shift indices with removal of annotations
                pam_start_idx = pam_start_idx - 1
                pam_end_idx = pam_end_idx - 2
                spacer_start_idx = 0
                spacer_end_idx = pam_start_idx

            else:
                pam_end_idx = pam_end_idx - 1
                cut_idx = cut_idx - 2
                spacer_start_idx = pam_end_idx
                spacer_end_idx = len(pe_format) - 3

        else:
            logger.error(
                'PE format parameter %s needs to cut site / within the spacer (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(
                    pe_format))
            sys.exit(1)

    else:
        logger.error(
            'PE format parameter %s needs to have one [PAM] present in its sequence (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(
                pe_format))
        sys.exit(1)

    # Remove annotations and convert into regex
    # 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
    # 'NNNNNNNNNNNNNNNNNNNNNGG'
    pe_format_rm_annotation = pe_format.replace('/', '').replace('[', '').replace(']', '')

    # Create PE format and PAM search sequences
    pe_format_search_plus = ''
    for base in pe_format_rm_annotation:
        pe_format_search_plus += iupac2bases(base)
    pe_format_search_minus = reverse_complement(pe_format_search_plus)

    pam_search = ''
    pam_sequence = pe_format_rm_annotation[pam_start_idx:pam_end_idx]

    for base in pam_sequence:
        pam_search += iupac2bases(base)

    ##### Initialize data storage for output
    pe_design = {}
    logger.info('Searching for pegRNAs for target sequences ...')
    counter = 1
    total_regions = len(target_design.keys())

    for target_name in target_design:
        # pegRNA spacer search for (+) and (-) strands with reference sequence
        reference_sequence = target_design[target_name]['reference_sequence']
        find_guides_ref_plus = [[m.start()] for m in
                                re.finditer('(?=%s)' % pe_format_search_plus, reference_sequence, re.IGNORECASE)]


        find_guides_ref_minus = [[m.start()] for m in
                         re.finditer('(?=%s)' % pe_format_search_minus, reference_sequence, re.IGNORECASE)]

        editnumber_sequence = target_design[target_name]['editnumber_sequence']

        find_guides_editnumber_plus = [[m.start()] for m in
                                       re.finditer('(?=%s)' % pam_search.replace('[', '[123456789'),
                                                   editnumber_sequence,
                                                   re.IGNORECASE)]

        find_guides_editnumber_minus = [[m.start()] for m in
                                        re.finditer(
                                            '(?=%s)' % reverse_complement(pam_search).replace('[', '[123456789'),
                                            editnumber_sequence, re.IGNORECASE)]

        editnumber2sequence = target_design[target_name]['editnumber2sequence']
        edit_sequence = target_design[target_name]['edit_sequence']
        edit_span_length = target_design[target_name]['edit_span_length'][1]
        is_shorter_than_30 = target_design[target_name]['is_shorter_than_30']


        # Find pegRNA spacers targeting (+) strand
        if find_guides_ref_plus:

            for match in find_guides_ref_plus:

                # Extract matched sequences and annotate type of prime editing
                full_search = reference_sequence[match[0]:match[0] + pe_format_length]
                spacer_sequence = full_search[spacer_start_idx:spacer_end_idx]
                extension_core_sequence = full_search[:cut_idx]
                downstream_sequence_ref = full_search[cut_idx:]
                downstream_sequence_length = len(downstream_sequence_ref)
                pam_ref = full_search[pam_start_idx:pam_end_idx]

                # Check to see if the extended non target strand is conserved in the edited strand
                try:
                    extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence,
                                                                                 edit_sequence).start(), re.search(
                        extension_core_sequence, edit_sequence).end()
                    downstream_sequence_edit = edit_sequence[
                                               extension_core_end_idx:extension_core_end_idx + downstream_sequence_length]
                    pam_edit = edit_sequence[extension_core_start_idx:extension_core_start_idx + pe_format_length][
                               pam_start_idx:pam_end_idx]

                    ## Annotate pegRNA
                    # Check if PAM is mutated relative to reference sequence
                    if pam_ref == pam_edit.upper():
                        pe_annotate = 'PAM_intact'

                    else:
                        # Check to see if mutation disrupts degenerate base positions within PAM
                        if re.search(pam_search, pam_edit.upper()):
                            pe_annotate = 'PAM_intact'

                        else:
                            pe_annotate = 'PAM_disrupted'

                    # Store pegRNA spacer
                    nick_ref_idx = match[0] + cut_idx
                    nick_edit_idx = extension_core_start_idx + cut_idx
                    target_design[target_name]['pegRNA']['+'].append(
                        [nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

                except:
                    continue


        # Find pegRNA spacers targeting (-) strand
        if find_guides_ref_minus:

            for match in find_guides_ref_minus:

                # Extract matched sequences and annotate type of prime editing
                full_search = reference_sequence[match[0]:match[0] + pe_format_length]
                spacer_sequence = full_search[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
                extension_core_sequence = full_search[pe_format_length - cut_idx:]
                downstream_sequence_ref = full_search[:pe_format_length - cut_idx]
                downstream_sequence_length = len(downstream_sequence_ref)
                pam_ref = full_search[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

                # Check to see if the extended non target strand is conserved in the edited strand
                try:
                    extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence,
                                                                                 edit_sequence).start(), re.search(
                        extension_core_sequence, edit_sequence).end()
                    downstream_sequence_edit = edit_sequence[
                                               extension_core_start_idx - downstream_sequence_length:extension_core_start_idx]
                    pam_edit = edit_sequence[extension_core_end_idx - pe_format_length:extension_core_end_idx][
                               pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

                    ## Annotate pegRNA
                    # Check if PAM is mutated relative to reference sequence
                    if pam_ref == pam_edit.upper():
                        pe_annotate = 'PAM_intact'

                    else:
                        # Check to see if mutation disrupts degenerate base positions within PAM
                        if re.search(reverse_complement(pam_search), pam_edit.upper()):
                            pe_annotate = 'PAM_intact'

                        else:
                            pe_annotate = 'PAM_disrupted'

                    # Store pegRNA spacer
                    nick_ref_idx = match[0] + (pe_format_length - cut_idx)
                    nick_edit_idx = extension_core_start_idx - downstream_sequence_length + (pe_format_length - cut_idx)
                    target_design[target_name]['pegRNA']['-'].append(
                        [nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

                except:
                    continue

        # Grab index information of edits to introduce to target sequence
        edit_start_in_ref = int(target_design[target_name]['edit_start_in_ref'])
        edit_stop_in_ref_rev = int(target_design[target_name]['edit_stop_in_ref_rev'])
        edit_span_length_w_ref = int(target_design[target_name]['edit_span_length'][0])
        edit_span_length_w_edit = int(target_design[target_name]['edit_span_length'][1])

        # Initialize pegRNA and ngRNA design dictionary
        pe_design[target_name] = {}

        ### regular pegRNA design here

        # Design pegRNAs targeting the (+) strand
        for peg_plus in target_design[target_name]['pegRNA']['+']:

            pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
            nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
            relative_index = nick2edit_length + 1
            pegid = '_'.join(
                map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', nick2edit_length]))
            pe_annotate_constant = pe_annotate

            # See if pegRNA spacer can introduce all edits
            if nick2edit_length >= 0 and editing_window_left[0] <= relative_index <= editing_window_left[1]:

                # Loop through RTT lengths
                for pbs_length in pbs_length_list:
                    pegRNA_ext = reverse_complement(
                        edit_sequence[
                        pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + edit_span_length_w_edit + nick2edit_length])
                    pegRNA_PBS = reverse_complement(
                        edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx])

                    pegRNA_RTT = reverse_complement(
                            edit_sequence[
                            pe_nick_edit_idx:pe_nick_edit_idx + edit_span_length_w_edit + nick2edit_length])

                    if pegid not in pe_design[target_name]:
                        # First list is for peg extension, second list is for nicking guide
                        pe_design[target_name][pegid] = [[]]

                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit


                    pe_design[target_name][pegid][0].append(
                            [pe_nick_ref_idx, pe_spacer_sequence,pe_pam_ref,
                             pe_annotate, '+', pbs_length,
                             pegRNA_ext, nick2lastedit_length, pegRNA_PBS, pegRNA_RTT, nick2edit_length,
                             edit_span_length_w_edit, pe_nick_edit_idx,pe_full_search])


        # Design pegRNAs targeting the (-) strand
        for peg_minus in target_design[target_name]['pegRNA']['-']:

            pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
            nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
            relative_index = nick2edit_length + 1
            pegid = '_'.join(
                map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-', nick2edit_length]))

            pe_annotate_constant = pe_annotate

            # See if pegRNA spacer can introduce all edits
            if nick2edit_length >= 0 and editing_window_right[0] <= relative_index <= editing_window_right[1]:


                for pbs_length in pbs_length_list:
                    pegRNA_ext = edit_sequence[
                                 pe_nick_edit_idx - edit_span_length_w_edit - nick2edit_length:pe_nick_edit_idx + pbs_length]
                    pegRNA_PBS = edit_sequence[pe_nick_edit_idx:pe_nick_edit_idx + pbs_length]
                    pegRNA_RTT = edit_sequence[
                                 pe_nick_edit_idx - edit_span_length_w_edit - nick2edit_length:pe_nick_edit_idx]

                    # Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
                    if pegid not in pe_design[target_name]:
                        pe_design[target_name][pegid] = [[]]

                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit



                    pe_design[target_name][pegid][0].append(
                                [pe_nick_ref_idx, reverse_complement(pe_spacer_sequence),reverse_complement(pe_pam_ref),
                                 pe_annotate, '-', pbs_length,
                                 pegRNA_ext, nick2lastedit_length, pegRNA_PBS, pegRNA_RTT, nick2edit_length,
                                 edit_span_length_w_edit, pe_nick_edit_idx,reverse_complement(pe_full_search)])


        if counter % 1000 == 0:
            logger.info('Completed pegRNA search for %s out of %s sites ...' % (counter, total_regions))
        counter += 1

    logger.info('Completed pegRNA search for %s out of %s sites ...' % (counter - 1, total_regions))

    # Output pegRNAs
    pegRNAs_summary_f = '%s_PrimeDesign.csv' % str(time.strftime("%Y%m%d_%I.%M.%S", time.localtime()))
    logger.info('Writing pegRNA designs into output file %s ...' % pegRNAs_summary_f)

    positive_candidates = []
    negative_candidates = []

    for target_name in pe_design:
        for pegid in pe_design[target_name]:

            strand = pegid.split('_')[-2]
            distance = int(pegid.split('_')[-1])

            if strand == '+':
                positive_candidates.append((pegid, distance))
            elif strand == '-':
                negative_candidates.append((pegid, distance))
    positive_candidates.sort(key=lambda x: x[1])
    negative_candidates.sort(key=lambda x: x[1])
    paired_designs = [(pos[0], neg[0]) for pos in positive_candidates for neg in negative_candidates]
    positive_dict = dict(positive_candidates)
    negative_dict = dict(negative_candidates)
    sorted_pairs = sorted(paired_designs, key=lambda x: (positive_dict[x[0]] + negative_dict[x[1]],
                                                         max(positive_dict[x[0]], negative_dict[x[1]])))


    def process_sequence_color(target_sequence, nick2edit_length_positive, nick2edit_length_negative, part2_length,
                         part1_length):
        # 步骤1：移除括号及其内部的字符
        left_bracket_index = target_sequence.find('(')
        right_bracket_index = target_sequence.find(')')  # 包括斜杠在内
        slash_index = target_sequence.find('/')
        left_sequence = target_sequence[:left_bracket_index]
        #print(target_sequence)
        #print(left_sequence)
        right_sequence = target_sequence[right_bracket_index + 1:]

        extracted_sequence = target_sequence[slash_index + 1:right_bracket_index]

        left_sequence_trimmed = left_sequence[
                                :-nick2edit_length_positive] if nick2edit_length_positive > 0 else left_sequence

        right_sequence_trimmed = right_sequence[
                                 nick2edit_length_negative:] if nick2edit_length_negative > 0 else right_sequence

        final_highlighted = left_sequence_trimmed + "<span style='color:red;'>" + extracted_sequence + "</span>" + right_sequence_trimmed



        return final_highlighted

    data_to_display = []
    counter = 1
    with open(out_dir + '/%s' % pegRNAs_summary_f, 'w') as f:
        # Headers for positive chain
        headers_positive = ['Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type',
                            'Spacer_sequence_positive', 'Spacer_GC_content_positive',
                            'PAM_sequence_positive', 'Extension_sequence_positive', 'Strand_positive',
                            'Annotation_positive', 'pegRNA-to-edit_distance_positive',
                            'Nick_index_positive', 'model_positive', 'PBS_length_positive',
                            'PBS_GC_content_positive', 'RTT_length_positive', 'RTT_GC_content_positive',
                            'First_extension_nucleotide_positive', 'PBS_positive', 'RTT_positive',
                            'PBS_Tm_positive', 'Recommend_Tm_positive', 'nick2edit_length_positive']

        # Headers for negative chain
        headers_negative = [
            'Spacer_sequence_negative', 'Spacer_GC_content_negative',
            'PAM_sequence_negative', 'Extension_sequence_negative', 'Strand_negative',
            'Annotation_negative', 'pegRNA-to-edit_distance_negative',
            'Nick_index_negative', 'PBS_length_negative',
            'PBS_GC_content_negative', 'RTT_length_negative', 'RTT_GC_content_negative',
            'First_extension_nucleotide_negative', 'PBS_negative', 'RTT_negative',
            'PBS_Tm_negative', 'Recommend_Tm_negative', 'nick2edit_length_negative','is_distance_issue']

        headers_combined = headers_positive + headers_negative
        headers_string = ','.join(map(str, headers_combined))
        f.write(headers_string + '\n')

        previous_Tms_positive = []
        previous_Tms_negative = []
        for target_name in pe_design:
            if sorted_pairs:
                for pair in sorted_pairs:
                    positive_pegid, negative_pegid = pair

                    for pegRNA_entry_positive in pe_design[target_name][positive_pegid][0]:
                        pe_nick_ref_idx_positive, pe_spacer_sequence_positive, pe_pam_ref_positive, pe_annotate_positive, pe_strand_positive, pbs_length_positive, pegRNA_ext_positive, nick2lastedit_length_positive, pegRNA_PBS_positive, pegRNA_RTT_positive, nick2edit_length_positive, edit_span_length_w_edit_positive, pe_nick_edit_idx_positive, full_search_positive = pegRNA_entry_positive

                        for pegRNA_entry_negative in pe_design[target_name][negative_pegid][0]:
                            pe_nick_ref_idx_negative, pe_spacer_sequence_negative, pe_pam_ref_negative, pe_annotate_negative, pe_strand_negative, pbs_length_negative, pegRNA_ext_negative, nick2lastedit_length_negative, pegRNA_PBS_negative, pegRNA_RTT_negative, nick2edit_length_negative, edit_span_length_w_edit_negative, pe_nick_edit_idx_negative, full_search_negative = pegRNA_entry_negative

                            assert edit_span_length_w_edit_positive == edit_span_length_w_edit_negative


                            is_distance_positive = nick2edit_length_positive > 0
                            is_distance_negative = nick2edit_length_negative > 0


                            is_distance_issue = is_distance_positive or is_distance_negative


                            total_length = nick2edit_length_positive + nick2edit_length_negative + edit_span_length_w_edit_positive


                            additional_length_needed = 0
                            is_shorter_than_30 = edit_span_length_w_edit_positive < homology_overlap


                            additional_length_needed = 0
                            if is_shorter_than_30:
                                additional_length_needed = homology_overlap - edit_span_length_w_edit_positive

                            if is_shorter_than_30:
                                avg_nick2edit_length = (nick2edit_length_positive + nick2edit_length_negative) // 2
                                part1_length = avg_nick2edit_length + additional_length_needed + edit_span_length_w_edit_positive
                                part2_length = avg_nick2edit_length + additional_length_needed + edit_span_length_w_edit_positive

                                # 如果 nick2edit 长度之和为奇数，将多出的1加到 part2_length 上
                                if (nick2edit_length_positive + nick2edit_length_negative) % 2 == 1:
                                    part2_length += 1

                            else:

                                base_part_length = (total_length - homology_overlap) // 2
                                part1_length = base_part_length + homology_overlap
                                part2_length = total_length - base_part_length

                                if part1_length + part2_length - homology_overlap != total_length:
                                    part2_length += 1

                            pbs_length_positive = len(pegRNA_PBS_positive)
                            pegRNA_ext_positive = reverse_complement(
                                edit_sequence[
                                pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length])
                            pegRNA_PBS_positive = reverse_complement(
                                edit_sequence[
                                pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive])

                            pegRNA_RTT_positive = reverse_complement(
                                edit_sequence[pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length])

                            if pegRNA_ext_positive[0].upper() == 'C' and pegRNA_ext_positive[0] is not None:

                                pegRNA_ext_positive = reverse_complement(
                                    edit_sequence[
                                    pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length + 1])
                                pegRNA_RTT_positive = reverse_complement(
                                    edit_sequence[pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length + 1])


                                if pegRNA_ext_positive[0].upper() == 'C' and pegRNA_ext_positive[0] is not None:
                                    pegRNA_ext_positive = reverse_complement(
                                        edit_sequence[
                                        pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length - 1])
                                    pegRNA_RTT_positive = reverse_complement(
                                        edit_sequence[
                                        pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length - 1])

                            pbs_length_negative = len(pegRNA_PBS_negative)
                            pegRNA_ext_negative = edit_sequence[
                                                  pe_nick_edit_idx_negative - part2_length:pe_nick_edit_idx_negative + pbs_length_negative]
                            pegRNA_PBS_negative = edit_sequence[
                                                  pe_nick_edit_idx_negative:pe_nick_edit_idx_negative + pbs_length_negative]
                            pegRNA_RTT_negative = edit_sequence[
                                                  pe_nick_edit_idx_negative - part2_length:pe_nick_edit_idx_negative]

                            if str(pegRNA_ext_negative[0]).upper() == 'C' and pegRNA_ext_negative[0] is not None:

                                pegRNA_ext_negative = edit_sequence[
                                                      pe_nick_edit_idx_negative - part2_length - 1:pe_nick_edit_idx_negative + pbs_length_negative]

                                pegRNA_RTT_negative = edit_sequence[
                                                      pe_nick_edit_idx_negative - part2_length - 1:pe_nick_edit_idx_negative]

                                if str(pegRNA_ext_negative[0]).upper() == 'C' and pegRNA_ext_negative[0] is not None:
                                    pegRNA_ext_negative = edit_sequence[
                                                          pe_nick_edit_idx_negative - part2_length + 1:pe_nick_edit_idx_negative + pbs_length_negative]

                                    pegRNA_RTT_negative = edit_sequence[
                                                          pe_nick_edit_idx_negative - part2_length + 1:pe_nick_edit_idx_negative]


                            pegRNA_ext_first_base_positive = pegRNA_ext_positive[0]
                            spacer_gc_content_positive = gc_content(pe_spacer_sequence_positive)
                            pbs_gc_content_positive = gc_content(pegRNA_PBS_positive)
                            rtt_gc_content_positive = gc_content(pegRNA_RTT_positive)


                            pegRNA_ext_first_base_negative = pegRNA_ext_negative[0]
                            spacer_gc_content_negative = gc_content(pe_spacer_sequence_negative)
                            pbs_gc_content_negative = gc_content(pegRNA_PBS_negative)
                            rtt_gc_content_negative = gc_content(pegRNA_RTT_negative)

                            if total_length % 2 == 1:
                                if pegRNA_ext_first_base_positive.upper() == 'C' or str(pegRNA_ext_first_base_negative).upper() == 'C':

                                    part1_length += 1
                                    part2_length -= 1
                                    pbs_length_positive = len(pegRNA_PBS_positive)
                                    pegRNA_ext_positive = reverse_complement(
                                        edit_sequence[
                                        pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length])
                                    pegRNA_PBS_positive = reverse_complement(
                                        edit_sequence[
                                        pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive])

                                    pegRNA_RTT_positive = reverse_complement(
                                        edit_sequence[
                                        pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length])

                                    if pegRNA_ext_positive[0].upper() == 'C':

                                        pegRNA_ext_positive = reverse_complement(
                                            edit_sequence[
                                            pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length + 1])
                                        pegRNA_RTT_positive = reverse_complement(
                                            edit_sequence[
                                            pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length + 1])

                                        if pegRNA_ext_positive[0].upper() == 'C':
                                            pegRNA_ext_positive = reverse_complement(
                                                edit_sequence[
                                                pe_nick_edit_idx_positive - pbs_length_positive:pe_nick_edit_idx_positive + part1_length - 1])
                                            pegRNA_RTT_positive = reverse_complement(
                                                edit_sequence[
                                                pe_nick_edit_idx_positive:pe_nick_edit_idx_positive + part1_length - 1])

                                    pbs_length_negative = len(pegRNA_PBS_negative)
                                    pegRNA_ext_negative = edit_sequence[
                                                          pe_nick_edit_idx_negative - part2_length:pe_nick_edit_idx_negative + pbs_length_negative]
                                    pegRNA_PBS_negative = edit_sequence[
                                                          pe_nick_edit_idx_negative:pe_nick_edit_idx_negative + pbs_length_negative]
                                    pegRNA_RTT_negative = edit_sequence[
                                                          pe_nick_edit_idx_negative - part2_length:pe_nick_edit_idx_negative]

                                    if str(pegRNA_ext_negative[0]).upper() == 'C':

                                        pegRNA_ext_negative = edit_sequence[
                                                              pe_nick_edit_idx_negative - part2_length - 1:pe_nick_edit_idx_negative + pbs_length_negative]

                                        pegRNA_RTT_negative = edit_sequence[
                                                              pe_nick_edit_idx_negative - part2_length - 1:pe_nick_edit_idx_negative]

                                        if str(pegRNA_ext_negative[0]).upper() == 'C':
                                            pegRNA_ext_negative = edit_sequence[
                                                                  pe_nick_edit_idx_negative - part2_length + 1:pe_nick_edit_idx_negative + pbs_length_negative]

                                            pegRNA_RTT_negative = edit_sequence[
                                                                  pe_nick_edit_idx_negative - part2_length + 1:pe_nick_edit_idx_negative]

                                    pegRNA_ext_first_base_positive = pegRNA_ext_positive[0]
                                    spacer_gc_content_positive = gc_content(pe_spacer_sequence_positive)
                                    pbs_gc_content_positive = gc_content(pegRNA_PBS_positive)
                                    rtt_gc_content_positive = gc_content(pegRNA_RTT_positive)


                                    pegRNA_ext_first_base_negative = pegRNA_ext_negative[0]
                                    spacer_gc_content_negative = gc_content(pe_spacer_sequence_negative)
                                    pbs_gc_content_negative = gc_content(pegRNA_PBS_negative)
                                    rtt_gc_content_negative = gc_content(pegRNA_RTT_negative)



                            pbs_Tm_positive = PBS_Tm(pegRNA_PBS_positive)
                            pbs_Tm_recommend_positive = recommend_Tm(pbs_Tm_positive,previous_Tms_positive)
                            previous_Tms_positive.append(pbs_Tm_positive)


                            pbs_Tm_negative = PBS_Tm(pegRNA_PBS_negative)
                            pbs_Tm_recommend_negative = recommend_Tm(pbs_Tm_negative,previous_Tms_negative)
                            previous_Tms_negative.append(pbs_Tm_negative)


                            final_seq = process_sequence_color(target_sequence, nick2edit_length_positive, nick2edit_length_negative, part2_length,part1_length)
                            pegRNA_RTT_positive = str(pegRNA_RTT_positive).upper()
                            pegRNA_RTT_negative = str(pegRNA_RTT_negative).upper()
                            pegRNA_ext_positive = str(pegRNA_ext_positive).upper()
                            pegRNA_ext_negative = str(pegRNA_ext_negative).upper()
                            if filter_tms:
                                if pbs_Tm_recommend_positive and pbs_Tm_recommend_negative:
                                    if Exclude_first_C:

                                        if not pegRNA_ext_first_base_positive.upper() == 'C' and not str(pegRNA_ext_first_base_negative).upper() == 'C':


                                            data_to_display.append([pe_spacer_sequence_positive, spacer_gc_content_positive,
                                                                    pe_pam_ref_positive,
                                                                    pegRNA_ext_positive, pe_strand_positive,
                                                                    pe_annotate_positive,
                                                                    nick2lastedit_length_positive,
                                                                    pe_nick_ref_idx_positive, 'TwinPE_model',
                                                                    pbs_length_positive,
                                                                    pbs_gc_content_positive,
                                                                    part1_length,
                                                                    rtt_gc_content_positive, pegRNA_ext_first_base_positive,
                                                                    pegRNA_PBS_positive,
                                                                    pegRNA_RTT_positive,
                                                                    pbs_Tm_positive, pbs_Tm_recommend_positive,
                                                                    nick2edit_length_positive,
                                                                    pe_spacer_sequence_negative, spacer_gc_content_negative,
                                                                    pe_pam_ref_negative,
                                                                    pegRNA_ext_negative, pe_strand_negative,
                                                                    pe_annotate_negative,
                                                                    nick2lastedit_length_negative,
                                                                    pe_nick_ref_idx_negative, pbs_length_negative,
                                                                    pbs_gc_content_negative,
                                                                    part2_length,
                                                                    rtt_gc_content_negative, pegRNA_ext_first_base_negative,
                                                                    pegRNA_PBS_negative,
                                                                    pegRNA_RTT_negative,
                                                                    pbs_Tm_negative, pbs_Tm_recommend_negative,
                                                                    nick2edit_length_negative, is_distance_issue, full_search_positive, full_search_negative,final_seq])
                                    else:

                                        data_to_display.append(
                                            [pe_spacer_sequence_positive, spacer_gc_content_positive,
                                             pe_pam_ref_positive,
                                             pegRNA_ext_positive, pe_strand_positive,
                                             pe_annotate_positive,
                                             nick2lastedit_length_positive,
                                             pe_nick_ref_idx_positive, 'TwinPE_model',
                                             pbs_length_positive,
                                             pbs_gc_content_positive,
                                             part1_length,
                                             rtt_gc_content_positive, pegRNA_ext_first_base_positive,
                                             pegRNA_PBS_positive,
                                             pegRNA_RTT_positive,
                                             pbs_Tm_positive, pbs_Tm_recommend_positive,
                                             nick2edit_length_positive,
                                             pe_spacer_sequence_negative, spacer_gc_content_negative,
                                             pe_pam_ref_negative,
                                             pegRNA_ext_negative, pe_strand_negative,
                                             pe_annotate_negative,
                                             nick2lastedit_length_negative,
                                             pe_nick_ref_idx_negative, pbs_length_negative,
                                             pbs_gc_content_negative,
                                             part2_length,
                                             rtt_gc_content_negative, pegRNA_ext_first_base_negative,
                                             pegRNA_PBS_negative,
                                             pegRNA_RTT_negative,
                                             pbs_Tm_negative, pbs_Tm_recommend_negative,
                                             nick2edit_length_negative, is_distance_issue, full_search_positive, full_search_negative,final_seq])


                            else:
                                if Exclude_first_C:
                                    if not pegRNA_ext_first_base_positive.upper() == 'C' and not str(pegRNA_ext_first_base_negative).upper() == 'C':

                                        data_to_display.append([pe_spacer_sequence_positive, spacer_gc_content_positive,
                                                                pe_pam_ref_positive,
                                                                pegRNA_ext_positive, pe_strand_positive,
                                                                pe_annotate_positive,
                                                                nick2lastedit_length_positive,
                                                                pe_nick_ref_idx_positive, 'TwinPE_model',
                                                                pbs_length_positive,
                                                                pbs_gc_content_positive,
                                                                part1_length,
                                                                rtt_gc_content_positive, pegRNA_ext_first_base_positive,
                                                                pegRNA_PBS_positive,
                                                                pegRNA_RTT_positive,
                                                                pbs_Tm_positive, pbs_Tm_recommend_positive,
                                                                nick2edit_length_positive,
                                                                pe_spacer_sequence_negative, spacer_gc_content_negative,
                                                                pe_pam_ref_negative,
                                                                pegRNA_ext_negative, pe_strand_negative,
                                                                pe_annotate_negative,
                                                                nick2lastedit_length_negative,
                                                                pe_nick_ref_idx_negative, pbs_length_negative,
                                                                pbs_gc_content_negative,
                                                                part2_length,
                                                                rtt_gc_content_negative, pegRNA_ext_first_base_negative,
                                                                pegRNA_PBS_negative,
                                                                pegRNA_RTT_negative,
                                                                pbs_Tm_negative, pbs_Tm_recommend_negative,
                                                                nick2edit_length_negative, is_distance_issue, full_search_positive, full_search_negative,final_seq])
                                else:

                                    data_to_display.append(
                                        [pe_spacer_sequence_positive, spacer_gc_content_positive,
                                         pe_pam_ref_positive,
                                         pegRNA_ext_positive, pe_strand_positive,
                                         pe_annotate_positive,
                                         nick2lastedit_length_positive,
                                         pe_nick_ref_idx_positive, 'TwinPE_model',
                                         pbs_length_positive,
                                         pbs_gc_content_positive,
                                         part1_length,
                                         rtt_gc_content_positive, pegRNA_ext_first_base_positive,
                                         pegRNA_PBS_positive,
                                         pegRNA_RTT_positive,
                                         pbs_Tm_positive, pbs_Tm_recommend_positive,
                                         nick2edit_length_positive,
                                         pe_spacer_sequence_negative, spacer_gc_content_negative,
                                         pe_pam_ref_negative,
                                         pegRNA_ext_negative, pe_strand_negative,
                                         pe_annotate_negative,
                                         nick2lastedit_length_negative,
                                         pe_nick_ref_idx_negative, pbs_length_negative,
                                         pbs_gc_content_negative,
                                         part2_length,
                                         rtt_gc_content_negative, pegRNA_ext_first_base_negative,
                                         pegRNA_PBS_negative,
                                         pegRNA_RTT_negative,
                                         pbs_Tm_negative, pbs_Tm_recommend_negative,
                                         nick2edit_length_negative, is_distance_issue, full_search_positive, full_search_negative,final_seq])

    return data_to_display









