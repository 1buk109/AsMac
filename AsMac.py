import argparse
import os
import time

import pandas as pd
import torch

from AsMac_model import AsMac
from AsMac_utility import *


def parse_args():
    parser = argparse.ArgumentParser(
        description='AsMac is a machine learning model for full length 16S rRNA sequence comparison')
    parser.add_argument("-i", "--input", required=True
                        , help="directory of the input fasta file)")
    parser.add_argument("-o", "--output", required=False
                        , help="prefered output directory. default: current working directory", default=os.getcwd())
    parser.add_argument("-m", "--model", required=False
                        , help="Model choice: 16S-full, 16S-V4, 16S-V3-V4, 23S-full, 23S-V5. default: full-16S", default='16S-full')

    return vars(parser.parse_args())


def parse_input_path(input_path):
    input_file = os.path.split(input_path)[1]
    input_file_name = input_file.split('.')[0]
    input_dir = os.path.split(input_path)[0]
    output_file_default = f'{input_file_name}_distances.csv'

    return input_file, input_file_name, input_dir, output_file_default


def parse_output_path(output_path):
    output_file = os.path.split(output_path)[1]
    output_file_name = output_file.split('.')[0]
    output_file = os.path.split(output_path)[1]
    output_file_name = output_file.split('.')[0]
    if output_path == os.getcwd():
        out_path = os.path.join(output_path, output_file_default)
    elif output_file.endswith(".csv"):
        out_path = output_path
    else:
        raise ValueError('Invalid output format, csv required.')

    return output_file, output_file_name, out_path


def read_fasta(seq_dir):
    seq_list = []
    info_list = []

    f_s = open(seq_dir, 'r')
    while 1:
        info = f_s.readline()[:-1]
        if not info:
            break
        seq = f_s.readline()[:-1]

        if set(seq).issubset({'A', 'T', 'U', 'G', 'C'}):
            info_list.append(info)
            seq_list.append(seq)
        else:
            raise ValueError('input sequence has symbol other than {A,T,U,C,G}')

    len_list = [len(s) for s in seq_list]
    max_l = max(len_list)
    min_l = min(len_list)

    return info_list, seq_list, max_l, min_l


def main():
    SEQ_LENGTHS = {'16S-full': 1400, '16S-V4': 151, '16S-V3-V4': 465, '23S-full': 3000, '23S-V5': 400}
    # define model
    EMBED_DIM = 300  # number of kernel sequences
    KERNEL_SIZE = 20  # kernel length

    args = parse_args()
    input_file, input_file_name, input_dir, output_file_default = parse_input_path(args["input"])
    output_file, output_file_name, out_path = parse_output_path(args["output"])
    model_type = os.path.split(args["model"])[1]

    if model_type in SEQ_LENGTHS.keys():
        print('Using model: %s, expecting sequence length of ~%i' % (model_type, SEQ_LENGTHS[model_type]))
    else:
        raise ValueError('Unsupported sequence type. Please choose one from : 16S-full, 16S-V4, 16S-V3-V4, 23S-full, 23S-V5.')

    # load sequences from input file
    info_list, seq_list, max_l, min_l = read_fasta(args["input"])
    print('%i sequences loaded. length range:[%i, %i]' % (len(seq_list), min_l, max_l))
    seq_oh = one_hot(seq_list)  # convert to one-hot

    # initialize AsMac model
    net = AsMac(4, EMBED_DIM, KERNEL_SIZE, model_type)
    net_state_dict = torch.load(f'model/{model_type}.pt')
    net.load_state_dict(net_state_dict)
    print('AsMac model loaded')

    # AsMac! Do the thing!
    tic = time.time()
    print('Computing embeddings for the sequences...')
    predictions = net.test_forward(seq_oh).detach().numpy().astype(np.float64)
    toc = time.time()
    print('Done!, cost %.2f seconds' % (toc - tic))
    np.fill_diagonal(predictions, 0)
    out_df = pd.DataFrame(predictions, columns=info_list, index=info_list)
    out_df.to_csv(out_path)
    print('Result saved in: %s' % out_path)


if __name__ == "__main__":
    main()
    