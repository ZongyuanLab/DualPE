from process import *
import argparse

def main():

    parser = argparse.ArgumentParser(description="Process gene editing parameters for Deletion.")

    # Editing type option (deletion, replacement, inversion)
    parser.add_argument('--file', type=str, required=True, help="File containing the target sequence")

    parser.add_argument('--edit_type', type=str, choices=['deletion', 'replacement', 'inversion'], default='deletion',
                        help="Type of edit (default: 'deletion')")

    parser.add_argument('--pam_sequence', type=str, default='NGG', help="PAM sequence (default: 'NGG')")

    parser.add_argument('--precise_option', action='store_true', help="Precise option (default: False)")
    parser.add_argument('--editing_window_left', type=int, default=15, help="Left editing window size (default: 15)")
    parser.add_argument('--editing_window_right', type=int, default=15, help="Right editing window size (default: 15)")
    parser.add_argument('--pbs_length_min', type=int, default=7, help="Minimum PBS length (default: 7)")
    parser.add_argument('--pbs_length_max', type=int, default=16, help="Maximum PBS length (default: 16)")
    parser.add_argument('--homology_overlap', type=int, default=30, help="Homology overlap (default: 30)")
    parser.add_argument('--pbs_Tm_Recommended', type=int, default=30, help="Recommended PBS Tm (default: 30)")
    parser.add_argument('--filter_tms', action='store_false', default=True, help="Filter TMS sequences (default: True)")

    parser.add_argument('--Exclude_first_C', action='store_false',default=True, help="Exclude first C (default: True)")


    # 解析参数
    args = parser.parse_args()
    with open(args.file, 'r') as f:
        target_sequence = f.read().strip()

    if args.edit_type == 'deletion':

        result = process_deletion_data(
            target_sequence = target_sequence,
            pam_sequence = args.pam_sequence,
            precise_option=args.precise_option,
            editing_window_left=args.editing_window_left,
            editing_window_right=args.editing_window_right,
            pbs_length_min=args.pbs_length_min,
            pbs_length_max=args.pbs_length_max,
            homology_overlap=args.homology_overlap,
            pbs_Tm_Recommended=args.pbs_Tm_Recommended,
            filter_tms=args.filter_tms,
            Exclude_first_C=args.Exclude_first_C,

        )

    elif args.edit_type == 'replacement':
        result = process_replacement_data(
            target_sequence=target_sequence,
            pam_sequence=args.pam_sequence,
            precise_option=args.precise_option,
            editing_window_left=args.editing_window_left,
            editing_window_right=args.editing_window_right,
            pbs_length_min=args.pbs_length_min,
            pbs_length_max=args.pbs_length_max,
            homology_overlap=args.homology_overlap,
            pbs_Tm_Recommended=args.pbs_Tm_Recommended,
            filter_tms=args.filter_tms,
            Exclude_first_C=args.Exclude_first_C,

        )

    elif args.edit_type == 'inversion':
        result = process_inversion_data(
            target_sequence=target_sequence,
            pam_sequence=args.pam_sequence,
            precise_option=args.precise_option,
            editing_window_left=args.editing_window_left,
            editing_window_right=args.editing_window_right,
            pbs_length_min=args.pbs_length_min,
            pbs_length_max=args.pbs_length_max,
            homology_overlap=args.homology_overlap,
            pbs_Tm_Recommended=args.pbs_Tm_Recommended,
            filter_tms=args.filter_tms,
            Exclude_first_C=args.Exclude_first_C,

        )







if __name__ == '__main__':
    main()