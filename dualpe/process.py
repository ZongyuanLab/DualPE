
import TwinPE_modelweb
import PD_modelweb
import invwebsim
from datetime import datetime

def adjust_sequence_with_unit(sequence, editing_window_left, editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find('/)')
   # print(unit_pos)
    # Adjust positions based on the editing window
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)
    new_unit_pos = max(0, unit_pos - editing_window_right - 1)  # +1 to treat '/)' as a unit
    #print(new_unit_pos)
    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace('(', '').replace('/)', '')

    # Re-insert the '(' at its new position

    sequence_adjusted = sequence_clean[:new_left_paren_pos] + '(' + sequence_clean[
                                                                        new_left_paren_pos:new_unit_pos] + '/)' + sequence_clean[
                                                                                                                  new_unit_pos:]

    return sequence_adjusted

def adjust_sequence_with_unit_1(sequence, editing_window_left):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find('/)')

    # Adjust positions based on the editing window
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)


    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace('(', '')


    sequence_adjusted = sequence_clean[:new_left_paren_pos] + '(' + sequence_clean[
                                                                        new_left_paren_pos:]

    return sequence_adjusted


def adjust_sequence_with_unit_2(sequence,editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find('/)')

    # Adjust positions based on the editing window

    new_unit_pos = max(0, unit_pos - editing_window_right)  # +1 to treat '/)' as a unit

    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace('/)', '')

    # Re-insert the '(' at its new position

    sequence_adjusted =  sequence_clean[:new_unit_pos] + '/)' + sequence_clean[new_unit_pos:]

    return sequence_adjusted

def generate_sequence_name():
    current_time = datetime.now()
    formatted_time = current_time.strftime("pesequence%Y%m%d%H%M%S")
    return formatted_time

def merge_unique_sublists_2(list1, list2, editing_window_left, editing_window_right):
    # 初始化集合和字典来跟踪唯一元素和最小值
    unique_keys_set = {(sublist[0], sublist[19]) for sublist in list1}
    unique_first_elements = {sublist[0] for sublist in list1}
    unique_twentieth_elements = {sublist[19] for sublist in list1}
    min_values_dict_18 = {sublist[0]: sublist[18] for sublist in list1 if len(sublist) > 19}
    min_values_dict_36 = {sublist[19]: sublist[36] for sublist in list1 if len(sublist) > 19}

    # 遍历第二个列表
    for sublist in list2:
        if len(sublist) > 19:  # 确保子列表长度足够
            key = (sublist[0], sublist[19])

            # 检查 sublist[0] 和 sublist[19] 是否之前未出现过
            if sublist[0] not in unique_first_elements:
                # sublist[18] = editing_window_left - sublist[18]  # 更新 sublist[18]
                unique_first_elements.add(sublist[0])  # 添加到集合中
                min_values_dict_18[sublist[0]] = sublist[18]  # 设置最小值

            if sublist[19] not in unique_twentieth_elements:
                # sublist[36] = editing_window_right - sublist[36]  # 更新 sublist[36]
                unique_twentieth_elements.add(sublist[19])  # 添加到集合中
                min_values_dict_36[sublist[19]] = sublist[36]  # 设置最小值

            # 检查子列表是否需要添加到list1
            if key not in unique_keys_set:
                # 更新或保留最小值
                if sublist[0] in min_values_dict_18 and sublist[18] > min_values_dict_18[sublist[0]]:
                    sublist[18] = min_values_dict_18[sublist[0]]
                if sublist[19] in min_values_dict_36 and sublist[36] > min_values_dict_36[sublist[19]]:
                    sublist[36] = min_values_dict_36[sublist[19]]

                # 将修改后的sublist添加到list1
                list1.append(sublist)
                # 更新集合以包含新元素
                unique_keys_set.add(key)

    return list1

def sort_sublists(list_to_sort):

    list_to_sort.sort(key=custom_sort)
    return list_to_sort

def custom_sort(sublist):

    sum_values = sublist[18] + sublist[36]

    max_value = max(sublist[18], sublist[36])

    return sum_values, -max_value

def select_corrected_deletion(editing_window_left,editing_window_right,sequence_data,sequence_data_gai,sequence_data_gai_1,sequence_data_gai_2):
    for data_gai in sequence_data_gai:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    for data_gai in sequence_data_gai_1:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

    for data_gai in sequence_data_gai_2:
        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    merged_list = merge_unique_sublists_2(sequence_data, sequence_data_gai, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_1, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_2, editing_window_left[1],
                                          editing_window_right[1])
    #print(merged_list[:5])
    sorted_list = filter_duplicate_data(merged_list)
    #print(sorted_list[:5])
    sorted_list = sort_sublists(sorted_list)

    return sorted_list

def filter_duplicate_data(sorted_list):

    data_dict = {}
    for data in sorted_list:
        key = (data[0], data[19])
        if key not in data_dict:
            data_dict[key] = []
        data_dict[key].append(data)


    new_sorted_list = []


    for key, datas in data_dict.items():
        if len(datas) > 1:

            to_keep = []


            for i in range(len(datas)):
                current = datas[i]
                keep_current = True

                for j in range(len(datas)):
                    if i != j:
                        other = datas[j]

                        if (current[16] + 2 == other[16]) or (current[34] + 2 == other[34]):
                            keep_current = False
                            break
                        elif (other[16] + 2 == current[16]) or (other[34] + 2 == current[34]):
                            keep_current = True
                            break

                if keep_current:
                    to_keep.append(current)


            new_sorted_list.extend(to_keep)
        else:

            new_sorted_list.extend(datas)

    return new_sorted_list

def get_pe_format(pam_sequence):
    if pam_sequence == 'NGG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGG]'

    elif pam_sequence == 'NRG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRG]'
    elif pam_sequence == 'SpG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGN]'
    elif pam_sequence == 'SpRYR':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRN]'
    elif pam_sequence == 'SpRYY':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NYN]'
    elif pam_sequence == 'NGCG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGCG]'
    elif pam_sequence == 'NGA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGA]'
    elif pam_sequence == 'NGT':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGT]'
    elif pam_sequence == 'NG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NG]'
    elif pam_sequence == 'NAAN':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NAAN]'
    elif pam_sequence == 'NNN':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNN]'
    elif pam_sequence == 'NRN':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRN]'
    elif pam_sequence == 'NGC':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGC]'
    elif pam_sequence == 'NRNH':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRNH]'
    elif pam_sequence == 'NGNG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGNG]'
    elif pam_sequence == 'NGN':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGN]'
    elif pam_sequence == 'NRCH':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRCH]'
    elif pam_sequence == 'NNGTGA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNGTGA]'
    elif pam_sequence == 'NNGG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNGG]'
    elif pam_sequence == 'NNAGAAW':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNAGAAW]'
    elif pam_sequence == 'NGGNG':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGGNG]'
    elif pam_sequence == 'NNNNGMTT':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNGMTT]'
    elif pam_sequence == 'NNNNCC':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNCC]'
    elif pam_sequence == 'NNGRRT':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNGRRT]'
    elif pam_sequence == 'NNNRRT':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNRRT]'
    elif pam_sequence == 'NNNVRYAC':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNVRYAC]'
    elif pam_sequence == 'NNNNRYAC':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNRYAC]'
    elif pam_sequence == 'NNNNGNA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNGNA]'
    elif pam_sequence == 'NRTA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRTA]'
    elif pam_sequence == 'NNNA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNA]'

    elif pam_sequence == 'NNNNCNR':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNCNR]'
    elif pam_sequence == 'NNNNCNAA':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NNNNCNAA]'
    elif pam_sequence == 'NRTH':
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NRTH]'

    elif pam_sequence == 'TBN':
        pe_format = '[TBN]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTTN':
        pe_format = '[TTTN]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTTV':
        pe_format = '[TTTV]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTNA':
        pe_format = '[TTNA]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTN':
        pe_format = '[TTN]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TG':
        pe_format = '[TG]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTCN':
        pe_format = '[TTCN]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TTTR':
        pe_format = '[TTTR]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'GGTT':
        pe_format = '[GGTT]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'KYTV':
        pe_format = '[KYTV]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'NTTR':
        pe_format = '[NTTR]/NNNNNNNNNNNNNNNNNNNNNNNN'


    elif pam_sequence == 'YTTV':
        pe_format = '[YTTV]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TYCV':
        pe_format = '[TYCV]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'TATV':
        pe_format = '[TATV]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'ATTN':
        pe_format = '[ATTN]/NNNNNNNNNNNNNNNNNNNNNNNN'
    elif pam_sequence == 'DTTN':
        pe_format = '[DTTN]/NNNNNNNNNNNNNNNNNNNNNNNN'

    return pe_format

def adjust_sequence_with_separated_slash_correct(sequence, editing_window_left, editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    slash_pos = sequence.find('/')
    right_paren_pos = sequence.find(')')

    # Adjust positions based on the editing windows
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)
    new_slash_pos = slash_pos - editing_window_right
    # new_right_paren_pos = right_paren_pos

    # Extract parts of the sequence
    # before_paren = sequence[:left_paren_pos]
    # after_paren_to_slash = sequence[left_paren_pos + 1:new_slash_pos]
    slash_to_paren = sequence[new_slash_pos:slash_pos]
    after_slash_to_right_paren = sequence[slash_pos + 1:right_paren_pos]
    after_right_paren = sequence[right_paren_pos + 1:]
    sequence = sequence.replace('(', '')
    # Construct the new sequence
    new_sequence = (
            sequence[:new_left_paren_pos] + '(' + sequence[new_left_paren_pos:new_slash_pos-1] + '/' + after_slash_to_right_paren + ')' +
            slash_to_paren + after_right_paren
    )

    return new_sequence

def adjust_sequence_with_separated_slash_correct_1(sequence, editing_window_left):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    slash_pos = sequence.find('/')
    right_paren_pos = sequence.find(')')

    # Adjust positions based on the editing windows
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)

    sequence = sequence.replace('(', '')
    # Construct the new sequence
    new_sequence = (
            sequence[:new_left_paren_pos] + '(' + sequence[new_left_paren_pos:]
    )

    return new_sequence

def adjust_sequence_with_separated_slash_correct_2(sequence, editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    slash_pos = sequence.find('/')
    right_paren_pos = sequence.find(')')

    # Adjust positions based on the editing windows

    new_slash_pos = slash_pos - editing_window_right

    slash_to_paren = sequence[new_slash_pos:slash_pos]
    after_slash_to_right_paren = sequence[slash_pos + 1:right_paren_pos]
    after_right_paren = sequence[right_paren_pos + 1:]

    # Construct the new sequence
    new_sequence = (sequence[:new_slash_pos] + '/' + after_slash_to_right_paren + ')' +
                    slash_to_paren + after_right_paren
                    )

    return new_sequence


def merge_unique_sublists_2(list1, list2, editing_window_left, editing_window_right):
    # 初始化集合和字典来跟踪唯一元素和最小值
    unique_keys_set = {(sublist[0], sublist[19]) for sublist in list1}
    unique_first_elements = {sublist[0] for sublist in list1}
    unique_twentieth_elements = {sublist[19] for sublist in list1}
    min_values_dict_18 = {sublist[0]: sublist[18] for sublist in list1 if len(sublist) > 19}
    min_values_dict_36 = {sublist[19]: sublist[36] for sublist in list1 if len(sublist) > 19}

    # 遍历第二个列表
    for sublist in list2:
        if len(sublist) > 19:
            key = (sublist[0], sublist[19])


            if sublist[0] not in unique_first_elements:
                # sublist[18] = editing_window_left - sublist[18]
                unique_first_elements.add(sublist[0])
                min_values_dict_18[sublist[0]] = sublist[18]

            if sublist[19] not in unique_twentieth_elements:

                unique_twentieth_elements.add(sublist[19])
                min_values_dict_36[sublist[19]] = sublist[36]


            if key not in unique_keys_set:

                if sublist[0] in min_values_dict_18 and sublist[18] > min_values_dict_18[sublist[0]]:
                    sublist[18] = min_values_dict_18[sublist[0]]
                if sublist[19] in min_values_dict_36 and sublist[36] > min_values_dict_36[sublist[19]]:
                    sublist[36] = min_values_dict_36[sublist[19]]


                list1.append(sublist)

                unique_keys_set.add(key)

    return list1
def select_corrected_replacement(editing_window_left,editing_window_right,sequence_data,sequence_data_gai,sequence_data_gai_1,sequence_data_gai_2):
    for data_gai in sequence_data_gai:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    for data_gai in sequence_data_gai_1:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

    for data_gai in sequence_data_gai_2:
        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    merged_list = merge_unique_sublists_2(sequence_data, sequence_data_gai, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_1, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_2, editing_window_left[1],
                                          editing_window_right[1])

    merged_list = filter_duplicate_data(merged_list)

    sorted_list = sort_sublists(merged_list)

    return sorted_list
def adjust_sequence_with_unit_inv(sequence, editing_window_left, editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find(')')
    print(unit_pos)
    # Adjust positions based on the editing window
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)
    new_unit_pos = max(0, unit_pos - editing_window_right - 1)  # +1 to treat '/)' as a unit
    print(new_unit_pos)
    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace('(', '').replace(')', '')

    # Re-insert the '(' at its new position

    sequence_adjusted = sequence_clean[:new_left_paren_pos] + '(' + sequence_clean[
                                                                    new_left_paren_pos:new_unit_pos] + ')' + sequence_clean[
                                                                                                             new_unit_pos:]

    return sequence_adjusted


def adjust_sequence_with_unit_1_inv(sequence, editing_window_left):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find('/)')

    # Adjust positions based on the editing window
    new_left_paren_pos = max(0, left_paren_pos + editing_window_left)

    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace('(', '')

    sequence_adjusted = sequence_clean[:new_left_paren_pos] + '(' + sequence_clean[
                                                                    new_left_paren_pos:]

    return sequence_adjusted


def adjust_sequence_with_unit_2_inv(sequence, editing_window_right):
    # Find the positions of special characters
    left_paren_pos = sequence.find('(')
    unit_pos = sequence.find(')')

    # Adjust positions based on the editing window

    new_unit_pos = max(0, unit_pos - editing_window_right)  # +1 to treat '/)' as a unit

    # Remove the original special characters from the sequence
    sequence_clean = sequence.replace(')', '')

    # Re-insert the '(' at its new position

    sequence_adjusted = sequence_clean[:new_unit_pos] + ')' + sequence_clean[new_unit_pos:]

    return sequence_adjusted


def select_corrected_inversion(editing_window_left, editing_window_right, sequence_data, sequence_data_gai,
                              sequence_data_gai_1, sequence_data_gai_2):
    for data_gai in sequence_data_gai:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    for data_gai in sequence_data_gai_1:
        if editing_window_left[1] < data_gai[18]:
            data_gai[18] = data_gai[18] - editing_window_left[1]
        else:
            data_gai[18] = editing_window_left[1] - data_gai[18]

    for data_gai in sequence_data_gai_2:
        if editing_window_right[1] < data_gai[36]:
            data_gai[36] = data_gai[36] - editing_window_right[1]
        else:
            data_gai[36] = editing_window_right[1] - data_gai[36]

    merged_list = merge_unique_sublists_2(sequence_data, sequence_data_gai, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_1, editing_window_left[1],
                                          editing_window_right[1])
    merged_list = merge_unique_sublists_2(merged_list, sequence_data_gai_2, editing_window_left[1],
                                          editing_window_right[1])
    sorted_list = []
    for data in merged_list:
        if data[18] <= editing_window_left[1] and data[36] <= editing_window_right[1]:
            sorted_list.append(data)
    sorted_list = filter_duplicate_data(sorted_list)
    sorted_list = sort_sublists(sorted_list)

    return sorted_list

def process_inversion_data(target_sequence,pam_sequence, precise_option=False, editing_window_left=15, editing_window_right=15,
                          pbs_length_min=7, pbs_length_max=16, homology_overlap=30, pbs_Tm_Recommended=30,
                          filter_tms=True, Exclude_first_C=True):

    editing_window_left = [1, editing_window_left]
    editing_window_right = [1, editing_window_right]
    pe_format = get_pe_format(pam_sequence)

    pbs_length_list = list(range(pbs_length_min, pbs_length_max + 1))
    target_name_primary = generate_sequence_name()

    sequence_data = invwebsim.main(target_name_primary, target_sequence, pe_format,
                                   editing_window_left, editing_window_right, pbs_length_list,
                                   homology_overlap,
                                   filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai = adjust_sequence_with_unit_inv(target_sequence, editing_window_left[1],
                                                        editing_window_right[1])
    sequence_data_gai = invwebsim.main(target_name_primary, target_sequence_gai, pe_format,
                                       editing_window_left, editing_window_right, pbs_length_list,
                                       homology_overlap,
                                       filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_1 = adjust_sequence_with_unit_1_inv(target_sequence, editing_window_left[1])
    sequence_data_gai_1 = invwebsim.main(target_name_primary, target_sequence_gai_1, pe_format,
                                         editing_window_left, editing_window_right, pbs_length_list,
                                         homology_overlap,
                                         filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_2 = adjust_sequence_with_unit_2_inv(target_sequence, editing_window_right[1])
    sequence_data_gai_2 = invwebsim.main(target_name_primary, target_sequence_gai_2, pe_format,
                                         editing_window_left, editing_window_right, pbs_length_list,
                                         homology_overlap,
                                         filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    sorted_list = select_corrected_inversion(editing_window_left, editing_window_right, sequence_data,
                                             sequence_data_gai,
                                             sequence_data_gai_1, sequence_data_gai_2)

    if not precise_option:
        if sorted_list:
            if sorted_list[0][18] == 0 and sorted_list[0][36] == 0:
                sorted_list = [sorted_list[0]]
            else:
                sorted_list = []
        else:
            sorted_list = []


    data_to_display = []
    for index, data in enumerate(sorted_list):
        sp_key = f"SP{index + 1}"



        result_dict = {
            'spacer_sequence_positive': data[0],
            'spacer_gc_content_positive': round(float(data[1]) * 100, 1),
            'pam_sequence_positive': data[2],
            'extension_sequence_positive': data[3],
            'strand_positive': data[4],
            'annotation_positive': data[5],
            'nick_to_edit_distance_positive': data[6],
            'nick_index_positive': data[7],
            'model_positive': data[8],
            'pbs_length_positive': data[9],
            'pbs_gc_content_positive': round(float(data[10]) * 100, 1),
            'rtt_length_positive': data[11],
            'rtt_gc_content_positive': round(float(data[12]) * 100, 1),
            'first_extension_nucleotide_positive': data[13],
            'pbs_positive': data[14],
            'rtt_positive': data[15],
            'pbs_tm_positive': data[16],
            'recommend_tm_positive': data[17],
            'nick2edit_length_positive': data[18],
            'spacer_sequence_negative': data[19],
            'spacer_gc_content_negative': round(float(data[20]) * 100, 1),
            'pam_sequence_negative': data[21],
            'extension_sequence_negative': data[22],
            'strand_negative': data[23],
            'annotation_negative': data[24],
            'nick_to_edit_distance_negative': data[25],
            'nick_index_negative': data[26],
            'pbs_length_negative': data[27],
            'pbs_gc_content_negative': round(float(data[28]) * 100, 1),
            'rtt_length_negative': data[29],
            'rtt_gc_content_negative': round(float(data[30]) * 100, 1),
            'first_extension_nucleotide_negative': data[31],
            'pbs_negative': data[32],
            'rtt_negative': data[33],
            'pbs_tm_negative': data[34],
            'recommend_tm_negative': data[35],
            'nick2edit_length_negative': data[36],


            'spacer_sequence_positive_length': len(data[0]),
            'spacer_sequence_negative_length': len(data[19]),
            'seq_before_inv': data[-2],
            'final_seq': data[-1],


        }
        data_to_display.append(result_dict)

    data_to_display = data_to_display[:20]
    display_results(data_to_display)

    return data_to_display

def process_replacement_data(target_sequence,pam_sequence, precise_option=False, editing_window_left=15, editing_window_right=15,
                          pbs_length_min=7, pbs_length_max=16, homology_overlap=30, pbs_Tm_Recommended=30,
                          filter_tms=True, Exclude_first_C=True):
    print(target_sequence)
    editing_window_left = [1, editing_window_left]
    editing_window_right = [1, editing_window_right]
    pe_format = get_pe_format(pam_sequence)

    pbs_length_list = list(range(pbs_length_min, pbs_length_max + 1))



    target_name_primary = generate_sequence_name()
    sequence_data = TwinPE_modelweb.main(target_name_primary, target_sequence, pe_format,
                                         editing_window_left, editing_window_right, pbs_length_list,
                                         homology_overlap,
                                         filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai = adjust_sequence_with_separated_slash_correct(target_sequence, editing_window_left[1],
                                                                       editing_window_right[1])
    sequence_data_gai = TwinPE_modelweb.main(target_name_primary, target_sequence_gai, pe_format,
                                             editing_window_left, editing_window_right, pbs_length_list,
                                             homology_overlap,
                                             filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_1 = adjust_sequence_with_separated_slash_correct_1(target_sequence, editing_window_left[1])
    sequence_data_gai_1 = TwinPE_modelweb.main(target_name_primary, target_sequence_gai_1, pe_format,
                                               editing_window_left, editing_window_right, pbs_length_list,
                                               homology_overlap,
                                               filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_2 = adjust_sequence_with_separated_slash_correct_2(target_sequence, editing_window_right[1])
    sequence_data_gai_2 = TwinPE_modelweb.main(target_name_primary, target_sequence_gai_2, pe_format,
                                               editing_window_left, editing_window_right, pbs_length_list,
                                               homology_overlap,
                                               filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    sorted_list = select_corrected_replacement(editing_window_left, editing_window_right, sequence_data,
                                               sequence_data_gai,
                                               sequence_data_gai_1, sequence_data_gai_2)
    if not precise_option:
        if sorted_list:
            if sorted_list[0][18] == 0 and sorted_list[0][36] == 0:
                sorted_list = [sorted_list[0]]
            else:
                sorted_list = []
        else:
            sorted_list = []
    data_to_display = []
    for index, data in enumerate(sorted_list):
        sp_key = f"SP{index + 1}"



        result_dict = {
            'spacer_sequence_positive': data[0],
            'spacer_gc_content_positive': round(float(data[1]) * 100, 1),
            'pam_sequence_positive': data[2],
            'extension_sequence_positive': data[3],
            'strand_positive': data[4],
            'annotation_positive': data[5],
            'nick_to_edit_distance_positive': data[6],
            'nick_index_positive': data[7],
            'model_positive': data[8],
            'pbs_length_positive': data[9],
            'pbs_gc_content_positive': round(float(data[10]) * 100, 1),
            'rtt_length_positive': data[11],
            'rtt_gc_content_positive': round(float(data[12]) * 100, 1),
            'first_extension_nucleotide_positive': data[13],
            'pbs_positive': data[14],
            'rtt_positive': data[15],
            'pbs_tm_positive': data[16],
            'recommend_tm_positive': data[17],
            'nick2edit_length_positive': data[18],
            'spacer_sequence_negative': data[19],
            'spacer_gc_content_negative': round(float(data[20]) * 100, 1),
            'pam_sequence_negative': data[21],
            'extension_sequence_negative': data[22],
            'strand_negative': data[23],
            'annotation_negative': data[24],
            'nick_to_edit_distance_negative': data[25],
            'nick_index_negative': data[26],
            'pbs_length_negative': data[27],
            'pbs_gc_content_negative': round(float(data[28]) * 100, 1),
            'rtt_length_negative': data[29],
            'rtt_gc_content_negative': round(float(data[30]) * 100, 1),
            'first_extension_nucleotide_negative': data[31],
            'pbs_negative': data[32],
            'rtt_negative': data[33],
            'pbs_tm_negative': data[34],
            'recommend_tm_negative': data[35],
            'nick2edit_length_negative': data[36],
            'is_distance_issue': data[37],

            'spacer_sequence_positive_length': len(data[0]),
            'spacer_sequence_negative_length': len(data[19]),
            'final_seq': data[-1],

        }

        data_to_display.append(result_dict)

    data_to_display = data_to_display[:20]
    display_results(data_to_display)

    return data_to_display

def process_deletion_data(target_sequence,pam_sequence, precise_option=False, editing_window_left=15, editing_window_right=15,
                          pbs_length_min=7, pbs_length_max=16, homology_overlap=30, pbs_Tm_Recommended=30,
                          filter_tms=True, Exclude_first_C=True):

    editing_window_left = [1, editing_window_left]
    editing_window_right = [1, editing_window_right]
    pe_format = get_pe_format(pam_sequence)

    pbs_length_list = list(range(pbs_length_min, pbs_length_max + 1))



    target_name_primary = generate_sequence_name()


    sequence_data = PD_modelweb.main(target_name_primary, target_sequence, pe_format,
                                     editing_window_left, editing_window_right, pbs_length_list,
                                     homology_overlap,
                                     filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai = adjust_sequence_with_unit(target_sequence, editing_window_left[1],
                                                    editing_window_right[1])
    sequence_data_gai = PD_modelweb.main(target_name_primary, target_sequence_gai, pe_format,
                                         editing_window_left, editing_window_right, pbs_length_list,
                                         homology_overlap,
                                         filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_1 = adjust_sequence_with_unit_1(target_sequence, editing_window_left[1])
    sequence_data_gai_1 = PD_modelweb.main(target_name_primary, target_sequence_gai_1, pe_format,
                                           editing_window_left, editing_window_right, pbs_length_list,
                                           homology_overlap,
                                           filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    target_sequence_gai_2 = adjust_sequence_with_unit_2(target_sequence, editing_window_right[1])
    sequence_data_gai_2 = PD_modelweb.main(target_name_primary, target_sequence_gai_2, pe_format,
                                           editing_window_left, editing_window_right, pbs_length_list,
                                           homology_overlap,
                                           filter_tms, Exclude_first_C, pbs_Tm_Recommended)

    sorted_list = select_corrected_deletion(editing_window_left, editing_window_right, sequence_data, sequence_data_gai,
                                            sequence_data_gai_1, sequence_data_gai_2)


    if not precise_option:
        if sorted_list:
            if sorted_list[0][18] == 0 and sorted_list[0][36] == 0:
                sorted_list = [sorted_list[0]]
            else:
                sorted_list = []
        else:
            sorted_list = []

    data_to_display = []
    for index, data in enumerate(sorted_list):
        sp_key = f"SP{index + 1}"

        result_dict = {
            'spacer_sequence_positive': data[0],
            'spacer_gc_content_positive': round(float(data[1]) * 100, 1),
            'pam_sequence_positive': data[2],
            'extension_sequence_positive': data[3],
            'strand_positive': data[4],
            'annotation_positive': data[5],
            'nick_to_edit_distance_positive': data[6],
            'nick_index_positive': data[7],
            'model_positive': data[8],
            'pbs_length_positive': data[9],
            'pbs_gc_content_positive': round(float(data[10]) * 100, 1),
            'rtt_length_positive': data[11],
            'rtt_gc_content_positive': round(float(data[12]) * 100, 1),
            'first_extension_nucleotide_positive': data[13],
            'pbs_positive': data[14],
            'rtt_positive': data[15],
            'pbs_tm_positive': data[16],
            'recommend_tm_positive': data[17],
            'nick2edit_length_positive': data[18],
            'spacer_sequence_negative': data[19],
            'spacer_gc_content_negative': round(float(data[20]) * 100, 1),
            'pam_sequence_negative': data[21],
            'extension_sequence_negative': data[22],
            'strand_negative': data[23],
            'annotation_negative': data[24],
            'nick_to_edit_distance_negative': data[25],
            'nick_index_negative': data[26],
            'pbs_length_negative': data[27],
            'pbs_gc_content_negative': round(float(data[28]) * 100, 1),
            'rtt_length_negative': data[29],
            'rtt_gc_content_negative': round(float(data[30]) * 100, 1),
            'first_extension_nucleotide_negative': data[31],
            'pbs_negative': data[32],
            'rtt_negative': data[33],
            'pbs_tm_negative': data[34],
            'recommend_tm_negative': data[35],
            'nick2edit_length_negative': data[36],


            'spacer_sequence_positive_length': len(data[0]),
            'spacer_sequence_negative_length': len(data[19]),
            'final_seq': data[-1],

        }

        data_to_display.append(result_dict)
    data_to_display = data_to_display[:20]
    display_results(data_to_display)

    return data_to_display


def display_results(data_to_display):
    #print(data_to_display)
    for index, result_dict in enumerate(data_to_display):
        print(f"Result {index + 1}:")

        print(
            f"  Spacer Sequence (Positive): {result_dict['spacer_sequence_positive']} | PAM Sequence: {result_dict['pam_sequence_positive']}")
        print(
            f"  PBS (Positive): {result_dict['pbs_positive']} | Length: {result_dict['pbs_length_positive']} | Tm: {result_dict['pbs_tm_positive']}℃")
        print(f"  RTT (Positive): {result_dict['rtt_positive']} | Length: {result_dict['rtt_length_positive']}")
        print(f"  Nick-to-Desired (Positive): {result_dict['nick2edit_length_positive']}")

        print(
            f"  Spacer Sequence (Negative): {result_dict['spacer_sequence_negative']} | PAM Sequence: {result_dict['pam_sequence_negative']}")
        print(
            f"  PBS (Negative): {result_dict['pbs_negative']} | Length: {result_dict['pbs_length_negative']} | Tm: {result_dict['pbs_tm_negative']}℃")
        print(f"  RTT (Negative): {result_dict['rtt_negative']} | Length: {result_dict['rtt_length_negative']}")
        print(f"  Nick-to-Desired (Negative): {result_dict['nick2edit_length_negative']}")

        print(f"  Sequence after edit: {result_dict['final_seq']}")

        print("\n" + "-" * 40 + "\n")







