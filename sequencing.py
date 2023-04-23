"""
Program: Clustering for True Induced Barcode Sequences
Author: Rameez Saiyid
Date: April 5, 2023

Description:
This script performs clustering for true induced barcode sequences. The input dataset
should be a list of RNA sequences in FASTQ format. The program will perform parallel 
clustering to group similar sequences together. The output will be a
text file containing the true sequences of the original dataset.

Usage:
The program can be run from the command line using the following command:
python sequencing.py
[filename] (10_5_candidates.fastq)
[cell count] (1570)
[anchor tag] (GTACTGCGGCCGCTACCTA)

True barcodes will be printed with confidence scores, and output to true_barcodes_output.txt
"""

import fastq as fq
import multiprocessing as mp
import matplotlib.pyplot as plt


def is_barcode(anchor):
    def is_it(f):
        return f.body[32:32+len(anchor)] == anchor
    return is_it


def get_barcode_prefix(seq):
    return seq.body[:30]


def group_barcodes_manual(barcodes):
    barcode_groups = []
    i = 0
    while barcodes:
        i += 1
        if i % 100000 == 0:
            print(i)
        barcode = barcodes.pop()
        # Find all the barcodes in the input list that differ from the current barcode by up to one character
        neighbors = [x for x in barcodes if sum(
            [a != b for a, b in zip(x, barcode)]) <= 1]
        # Add the current barcode and its neighbors to a new group
        barcode_groups.append(neighbors + [barcode])
        # Remove the barcodes in the current group from the input list
        pregroup_count = (len(barcodes))
        barcodes = [barcode for barcode in barcodes if barcode not in neighbors]
        postgroup_count = (len(barcodes))
        print("grouped", pregroup_count - postgroup_count)

    return barcode_groups


def find_neighbors(barcode, barcodes):
    return [x for x in barcodes if sum([a != b for a, b in zip(x, barcode)]) <= 1]


def group_barcodes_parallel(barcodes):
    barcode_groups = []
    barcode_set = set(barcodes)
    processed_barcodes = set()
    pool = mp.Pool()
    last_status = 0

    while barcode_set:
        barcode = barcode_set.pop()
        # Skip if barcode is already processed
        if barcode in processed_barcodes:
            continue

        # Use multiprocessing to find neighbors
        neighbors = pool.starmap(find_neighbors, [(barcode, barcode_set)])
        # Flatten list of neighbors
        neighbors = [item for sublist in neighbors for item in sublist]
        barcode_group = neighbors + [barcode]
        # Add the barcode group to the list of groups
        barcode_groups.append(barcode_group)
        # Update the set of processed barcodes
        processed_barcodes.update(barcode_group)

        # Remove processed barcodes from the barcode set
        barcode_set.difference_update(barcode_group)
        if (cur_percent := (float(len(processed_barcodes)) / float(len(barcode_set) + len(processed_barcodes)))*100) - last_status > 10:
            last_status = cur_percent
            print(f"{format(cur_percent, '.2f')}% clustered")

    print(f"100% clustered")
    pool.close()

    return barcode_groups


def generate_bar_chart(x, y, x_title, y_title, title):
    # create a bar chart
    plt.bar(x, y)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.title(title)
    plt.show()


if __name__ == '__main__':
    file_name = input("Please enter fastq filename: ")
    num_cells = input(
        "How many cells were sequenced? (Press enter to default to 100) ")
    if not len(num_cells):
        num_cells = 100
    anchor = input("What is the known anchor sequence? ")
    if not len(anchor):
        # Default if none provided
        anchor = "GTACTGCGGCCGCTACCTA"
    processing_type = "fast"
    # 10_6_candidates.fastq

    raw_data = []
    try:
        fos = fq.read(file_name)
        raw_data = list(fos)
    except Exception as e:
        print(e)

    potentials = list(filter(is_barcode(anchor), raw_data))
    barcode_potentials = list(map(get_barcode_prefix, potentials))
    print(f"Analyzing {len(barcode_potentials)} potential barcodes")

    # parallel clustering:
    grouped_barcodes = []
    if processing_type == "fast":
        grouped_barcodes = group_barcodes_parallel(
            barcode_potentials
        )
    else:
        # Takes 8 hours...
        grouped_barcodes = group_barcodes_manual(barcode_potentials)

    total_count = sum(list(map(len, grouped_barcodes)))
    print(f"There are {total_count} potentially true barcodes")

    # Process group and variance counts
    # Sort groups by decreasing size (higher confidence true)
    true_barcodes = []
    barcode_details = []
    for group in range(len(grouped_barcodes)):
        barcode_group = []
        for barcode in grouped_barcodes[group]:
            barcode_group.append({
                "barcode": barcode,
                "count": barcode_potentials.count(barcode)
            })
        barcode_group.sort(key=lambda variant: variant["count"], reverse=True)
        barcode_details.append(barcode_group)

    barcode_details.sort(key=lambda group: sum(
        list(map(lambda seq: seq["count"], group))), reverse=True)

    # =========ANALYTICS=========
    # group_sums = {group: sum(
    #     list(map(lambda seq: seq["count"], barcode_details[group]))) for group in range(len(barcode_details))}
    # generate_bar_chart(group_sums.keys(), group_sums.values(
    # ), 'Groups of Sequences', 'Occurances', 'Sequence Counts Across Groups')

    # first_group_sums = {variant["barcode"]: variant["count"]
    #                     for variant in barcode_details[0]}
    # generate_bar_chart(first_group_sums.keys(), first_group_sums.values(
    # ), 'Variants', 'Occurances', 'Variance Counts for a Sequence')
    # ===========================
    # ===========================

    for group in range(min(int(num_cells), len(barcode_details))):
        barcode_group = barcode_details[group]
        most_freq_seq = max(barcode_group, key=lambda seq: seq["count"])
        total_group_count = sum(list(map(lambda x: x["count"], barcode_group)))
        # Calculate confidence of barcode variant amongst group of hamming-1 neighbors by frequency
        most_freq_seq["confidence"] = format(
            float(most_freq_seq["count"]) / float(total_group_count), '.2f')
        true_barcodes.append(most_freq_seq)
    true_barcodes.sort(key=lambda x: x["confidence"], reverse=True)

    print(
        f"=================== {len(true_barcodes)} True Barcodes ===================")
    with open("true_barcodes_output.txt", "w") as file:
        for barcode in true_barcodes:
            print(barcode)
            file.write(barcode["barcode"] + "\n")
