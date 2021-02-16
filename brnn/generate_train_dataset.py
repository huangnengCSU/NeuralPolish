"""
 * @Author: huangneng 
 * @Date: 2019-12-05 15:26:11 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2019-12-05 15:26:11 
 """

import os
import subprocess
import pysam
import argparse
import shutil
from parse_mPileup_to_feature import parse_mPileup
from multiprocessing import Pool
from itertools import repeat
import time


def get_mPileup_regions(sam_file, win_size=64):
    sam = pysam.AlignmentFile(sam_file)
    n_refernece = sam.nreferences
    references = sam.references
    ref_lengths = sam.lengths

    mpileup_region_list = []
    for i in range(n_refernece):
        contig = references[i]
        ctg_len = ref_lengths[i]
        start_pos = 1
        """
        In samtools mpileup output file, contig start with 1, instead of 0
        input region (chromsome:p1-p2), both p1 and p2 position are contained \
            in the output file
        """
        while start_pos + (win_size - 1) <= ctg_len:
            end_pos = start_pos + (win_size - 1)  # start=1, end=64, win_size=64
            region_str = "%s:%s-%s" % (
                str(contig),
                str(start_pos),
                str(end_pos),
            )  # chromosome:1-64
            mpileup_region_list.append(region_str)
            start_pos = start_pos + win_size  # new_start=65
        if start_pos <= ctg_len:
            end_pos = ctg_len
            region_str = "%s:%s-%s" % (str(contig), str(start_pos), str(end_pos))
            mpileup_region_list.append(region_str)
    return mpileup_region_list


def run_samtools_and_parse_alignment(
    ref_file, bam_file, region, output_dir, win_size, max_insertion_size, mode
):
    mpileup_output = output_dir + "/" + region + "." + mode
    samtools_command = "samtools mpileup -B --reference {reference} -r {region_str} --output-QNAME {bamfile} -o {output_file} >> log.txt 2>&1".format(
        reference=ref_file,
        region_str=region,
        bamfile=bam_file,
        output_file=mpileup_output,
    )
    # print(samtools_command)
    returnCode = subprocess.run(samtools_command, shell=True, check=True)
    if returnCode.returncode==0:
        parse_mPileup(mpileup_output, output_dir, win_size, max_insertion_size, mode)
        os.remove(mpileup_output)
    else:
        print("Program execution encountered an error! exit.")
        os._exit(1)


def threads_generate_training_features(
    ref_file,
    bam_file,
    feature_output_dir,
    threads=2,
    win_size=64,
    max_insertion_size=5,
    mode="feature",
):
    region_list = get_mPileup_regions(bam_file, win_size)
    pool = Pool(threads)
    thread_params = zip(
        repeat(ref_file),
        repeat(bam_file),
        region_list,
        repeat(feature_output_dir),
        repeat(win_size),
        repeat(max_insertion_size),
        repeat(mode),
    )
    pool.starmap(run_samtools_and_parse_alignment, thread_params)


if __name__ == "__main__":
    start_time = time.time()
    parse = argparse.ArgumentParser()
    parse.add_argument("-r", "--ref", help="reference file", required=True)
    parse.add_argument("-b", "--bam", help="alignment bam file", required=True)
    parse.add_argument(
        "-o",
        "--output",
        help="directory of feature output or label output",
        required=True,
    )
    parse.add_argument("-t", "--threads", help="number of threads", default=4, type=int)
    parse.add_argument(
        "-w", "--window_size", help="window size of each sample", default=64, type=int
    )
    parse.add_argument(
        "-e",
        "--max_insert_size",
        help="max insertion size to be encoded",
        default=5,
        type=int,
    )
    parse.add_argument(
        "-m",
        "--mode",
        help="feature mode will reserve reference feature, label mode don't",
        choices=["feature", "label"],
        default="feature",
    )

    args = parse.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    else:
        rm_confirm = "no"
        rm_confirm = input(
            "do you want to delete the directory %s \n please input `yes` or `no`:"
            % args.output
        )
        if rm_confirm == "yes":
            shutil.rmtree(args.output)
            print("remove directory %s completed." % args.output)
            os.mkdir(args.output)

    threads_generate_training_features(
        args.ref,
        args.bam,
        args.output,
        args.threads,
        args.window_size,
        args.max_insert_size,
        args.mode,
    )
    print(
        "generate_train_dataset.py total time cost: %.5f sec"
        % (time.time() - start_time)
    )

