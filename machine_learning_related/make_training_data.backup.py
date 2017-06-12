"""
Convert the raw sequence and the lables to npy arrays for faster batch reading.
Split data into training, test and validation set.
Will store a .npz file with the labels and sequences and a coord file per test/valid and train set
"""
# from __future__ import absolute_import
# from __future__ import division
# from __future__ import print_function

import numpy as np
import h5py
import argparse
from operator import itemgetter

# Define arguments -------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="""Take a raw sequence and a labels bed like file and encode and store
    both as numpy arrays. Split up into traiing, test and validation samples.""")
parser.add_argument('in_file', type=str,
                    help='Five column file. [chr start end comma separated IDs to split and raw sequence].')
parser.add_argument('--split_mode', dest='split_mode', default='random', choices=['random', 'chr'],
    help="""Specify how to split up the data into training, test
    and validation set. []chr] - select chromomes from which
    features are attritbuted to the different sets. [random]
    - split the features by random sampling. Needs --chr_test
    --chr_valid or --frac_test --frac_valid declared respectively.
    Default = [random]""")
parser.add_argument('--frac_test', type=float, dest='frac_test', default=0.05,
    help='Fraction of total set to sample into test set. (Float > 1.0)')
parser.add_argument('--frac_valid', type=float, dest='frac_valid', default=0.05,
    help='Fraction of total set to sample into validation set. (Float > 1.0)')
parser.add_argument('--chr_test', nargs='+', dest='chr_test', default='chr20',
    help="""Select one ore more space separated chromosome (chr1 chr2 chr3) to use
     as test chromosomes. Default = chr20 Only if split_mode = 'chr' """)
parser.add_argument('--chr_valid', nargs='+', dest='chr_valid', default='chr21',
    help="""Select one ore more space separated chromosome (chr1 chr2 chr3) to use
     as validation chromosomes. Default = chr20 Only if split_mode = 'chr' """)
parser.add_argument('--save_prefix', dest='save_prefix', default='./data_set',
    help='Prefix to store the training/ test and validation sets. Default = ./data_set')
parser.add_argument('--seed', dest='seed', type=int, default=1234,
    help='Random seed for sampling.')
# Parse arguments
args = parser.parse_args()

# Helper get hotcoded sequence
def get_hot_coded_seq(sequence):
    """Convert a 4 base letter sequence to 4-row x-cols hot coded sequence"""
    # initialise empty
    hotsequence = np.zeros((len(sequence),4))
    # set hot code 1 according to gathered sequence
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            hotsequence[i,0] = 1
        elif sequence[i] == 'C':
            hotsequence[i,1] = 1
        elif sequence[i] == 'G':
            hotsequence[i,2] = 1
        elif sequence[i] == 'T':
            hotsequence[i,3] = 1
    # return the numpy array
    return hotsequence

print("\n# === Creating a Training, Test and Validation Set from provided input === #")

# Sed seed for random sampling -------------------------------------------------
np.random.seed(args.seed)

# Read in data -----------------------------------------------------------------
print("\nReading lines ...")

# read data in, split into vectors
lines = open(args.in_file, "r").readlines()
chroms = []
start = []
stop = []
label = []
label_tmp = []
seq = []
for l in lines:
    l = l.rstrip()
    l = l.split("\t")
    chroms.append(l[0])
    start.append(l[1])
    stop.append(l[2])
    label.append(l[3])
    label_tmp.append(l[3].split(","))
    seq.append(l[4])

# Get all lables and sort and assign them to binary labels ---------------------
label_tmp = [item for sublist in label_tmp for item in sublist]
label_tmp = np.array(label_tmp)
unique_labels = np.unique(label_tmp)
print(unique_labels)
num_ids = len(unique_labels)  # get number of unique ids
print("\nNumber of distinct labels found: " + str(num_ids))
print("\nDistinct labels: " + ' '.join(map(str,unique_labels)))
# init binary representatons --------------------------------------------------
# make a look-up dictionary with a binary label per id
bin_look_up = {}
for i in range(num_ids):
    bin_look_up[unique_labels[i]] = np.zeros((num_ids), dtype=np.uint64)
    bin_look_up[unique_labels[i]][i] = 1
# print a table with the intial labels for future reference
print("\nConverting to binary representation:")
for i in range(num_ids):
    print(unique_labels[i] + " -->\t" + ','.join(map(str, bin_look_up[unique_labels[i]])))
# Go through labels per seq and sum up a binary representing all active IDs ----
label_bin = np.zeros((len(label), num_ids),  dtype=np.float)
for j in range(len(label)):
    l = label[j].split(",") # split by commat
    for i in range(len(l)):
        label_bin[j,] = label_bin[j,] + bin_look_up[l[i]]

# Sample Test/ Validation and Training set according to selected mode -----------
input_rows = np.array(range(len(chroms)))  # make an array of input rows to sample from once

# if to split based on fractions randomly form all chromosomes
if args.split_mode == 'random':
    print("\nSampling randomly across chromosomes.")
    to_sample_test = round(len(input_rows) * args.frac_test)  # get fractions
    to_sample_valid = round(len(input_rows) * args.frac_valid)
    to_sample_train = len(input_rows) - to_sample_test - to_sample_valid
    print("%s Test cases\n%s Validation cases\n%s Training cases left." %
        (int(to_sample_test), int(to_sample_valid), int(to_sample_train)))
    # sample and get test and valid rows
    tmp_sampled = np.random.choice(input_rows, size=int(to_sample_test+to_sample_valid), replace=False)
    test_rows = tmp_sampled[range(int(to_sample_test))]
    valid_rows = tmp_sampled[range(int(to_sample_test), int(to_sample_test+to_sample_valid))]
    # prune remaining training cases
    training_rows = np.delete(input_rows, tmp_sampled)
    # random resample training rows
    training_rows = np.random.choice(training_rows, training_rows.size, replace=False)


elif args.split_mode == 'chr':
    print("\nSetting specifc chromosomes as test and validatipon set:")
    print("Using %s as Test, %s as Validation and the remaining as Training cases" % (args.chr_test, args.chr_valid))

    test_rows = []
    valid_rows = []

    # match row numbers against chromosomes
    for i in range(len(chroms)):
        if chroms[i] in args.chr_test:
            test_rows.append(i)
        if chroms[i] in args.chr_valid:
            valid_rows.append(i)
    # prune remaining training cases
    training_rows = np.delete(input_rows, (test_rows + valid_rows))
    print("%s Test cases \n%s Validation cases\n%s Training cases left." %
        (int(len(test_rows)), int(len(valid_rows)), int(len(training_rows))))
    # Resample for randomness even in chromosome case
    training_rows = np.random.choice(training_rows, training_rows.size, replace=False)
    test_rows = np.random.choice(test_rows, test_rows.size, replace=False)
    valid_rows = np.random.choice(valid_rows, valid_rows.size, replace=False)

# Convert sequences to hot coded numpy array
print("\nConverting sequences to hot coded representations ...")
seq_hot = np.empty((len(seq), len(seq[0]), 4), dtype=np.float)
for i in range(len(seq)):
    seq_hot[i,:] = get_hot_coded_seq(seq[i])
del seq
print(seq_hot.shape)


# TODO read without converting the sequence -- determine training and test rows etc than read the input file again and convert and assign the labels and sequence and directly write to  output file without storing everything


print("\nStoring ...")
# write training/test/validation set coords ------------------------------------
write_train_coords = open(args.save_prefix + "_training_coords.bed", "w")
for tr in training_rows:
    write_train_coords.write("%s\t%s\t%s\n" % (chroms[tr], start[tr], stop[tr]))
write_test_coords = open(args.save_prefix + "_test_coords.bed", "w")
for tr in test_rows:
    write_test_coords.write("%s\t%s\t%s\n" % (chroms[tr], start[tr], stop[tr]))
write_valid_coords = open(args.save_prefix + "_validation_coords.bed", "w")
for tr in valid_rows:
    write_valid_coords.write("%s\t%s\t%s\n" % (chroms[tr], start[tr], stop[tr]))

# Save Labels and Sequences in comboned npz file
# Save Training Labels and Sequences
training_labels = label_bin[training_rows,]
training_seqs = seq_hot[training_rows,]
# Save Test Labels and Sequences
test_labels = label_bin[test_rows,]
test_seqs = seq_hot[test_rows,]

# save traiing data in hdf5
h5f = h5py.File(args.save_prefix + "_training_data.h5", 'w')
h5f.create_dataset('training_seqs', data=training_seqs)
h5f.create_dataset('training_labels', data=training_labels)
h5f.create_dataset('test_seqs', data=test_seqs)
h5f.create_dataset('test_labels', data=test_labels)
h5f.close()

# Save Validation Labels and Sequences
validation_labels = label_bin[valid_rows,]
validation_seqs = seq_hot[valid_rows,]
h5f = h5py.File(args.save_prefix + "_validation_data.h5", 'w')
h5f.create_dataset('validation_seqs', data=validation_seqs)
h5f.create_dataset('validation_labels', data=validation_labels)
h5f.close()

print("\nSuccessfully saved a coord and a *data.npz file for each subset.\n")
