"""
 * @Author: huangneng 
 * @Date: 2019-12-19 09:53:48 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2019-12-19 09:53:48 
"""

import os
import re
import sys
import argparse
import time


class alignment_record_sam:
    Qname = ""
    flag = 0
    Rname = ""
    rs = -1
    map_quality = 0
    cigar = 0

    def get_mapping_length(self, max_indels):
        pattern = re.compile(r"((\d)+(X|=|I|D))")
        it = pattern.finditer(self.cigar)
        ref_len = 0
        for match in it:
            if match.group().endswith("="):
                length = int(match.group(0)[:len(match.group(0)) - 1])
                ref_len += length
            elif match.group().endswith("X"):
                length = int(match.group(0)[:len(match.group(0)) - 1])
                ref_len += length
            elif match.group().endswith("D"):
                length = int(match.group(0)[:len(match.group(0)) - 1])
                if length > max_indels:
                    ref_len = -1
                    break
                ref_len += length
            elif match.group().endswith("I"):
                length = int(match.group(0)[:len(match.group(0)) - 1])
                if length > max_indels:
                    ref_len = -1
                    break
            else:
                continue
        return ref_len


def get_alignment_record_with_bar(sam, max_indels, min_baseq, line_num):
    record_dict = {}
    fsam = open(sam, "r")
    line_idx = 0
    for line in fsam:
        if line.startswith("@"):
            line_idx += 1
            continue
        else:
            line_idx += 1
            records = alignment_record_sam()
            line = line.strip("\n")
            line_items = line.split("\t")
            records.Qname = line_items[0]
            records.flag = int(line_items[1])
            records.Rname = line_items[2]
            records.rs = int(line_items[3]) - 1
            records.map_quality = int(line_items[4])
            records.cigar = line_items[5]
            if records.map_quality < min_baseq:
                continue
            mapping_length = records.get_mapping_length(max_indels)
            record_dict[(records.Qname, records.Rname)] = record_dict.get(
                (records.Qname, records.Rname), [-1, -1])  # 比对长度，行号
            if mapping_length > record_dict[(records.Qname, records.Rname)][0]:
                record_dict[(records.Qname, records.Rname)] = [
                    mapping_length,
                    line_idx,
                ]  # 更新
            print("\r process %d / %d" % (line_idx, line_num), end="")
    fsam.close()
    print("")
    print("original sam file records number: %d" % line_idx)
    return record_dict


def get_alignment_record(sam, max_indels, min_baseq):
    record_dict = {}
    fsam = open(sam, "r")
    line_idx = 0
    for line in fsam:
        if line.startswith("@"):
            line_idx += 1
            continue
        else:
            line_idx += 1
            records = alignment_record_sam()
            line = line.strip("\n")
            line_items = line.split("\t")
            records.Qname = line_items[0]
            records.flag = int(line_items[1])
            records.Rname = line_items[2]
            records.rs = int(line_items[3]) - 1
            records.map_quality = int(line_items[4])
            records.cigar = line_items[5]
            if records.map_quality < min_baseq:
                continue
            mapping_length = records.get_mapping_length(max_indels)
            record_dict[(records.Qname, records.Rname)] = record_dict.get(
                (records.Qname, records.Rname), [-1, -1])  # 比对长度，行号
            # print(records.Qname, records.Rname, mapping_length)
            if mapping_length > record_dict[(records.Qname, records.Rname)][0]:
                record_dict[(records.Qname, records.Rname)] = [
                    mapping_length,
                    line_idx,
                ]  # 更新
    fsam.close()
    print("original sam file records number: %d" % line_idx)
    return record_dict


def filter_sam(sam, outfile, max_indels, min_baseq, line_num):
    if line_num is not None:
        Progress_bar = True
    else:
        Progress_bar = False
    if Progress_bar:
        record_dict = get_alignment_record_with_bar(sam, max_indels, min_baseq,
                                                    line_num)
    else:
        record_dict = get_alignment_record(sam, max_indels, min_baseq)
    reserve_line_idxs = []
    for key in record_dict.keys():
        reserve_line_idxs.append(record_dict[key][1])
    print("reserved file records number: %d" % len(reserve_line_idxs))
    # print(reserve_line_idxs)
    reserve_line_idxs = set(reserve_line_idxs)
    fout = open(outfile, "w")
    with open(sam) as fsam:
        line_idx = 0
        for line in fsam:
            if line.startswith("@"):
                line_idx += 1
                fout.write(line)
            else:
                line_idx += 1
                if line_idx in reserve_line_idxs:
                    fout.write(line)
    fout.close()


if __name__ == "__main__":
    start_time = time.time()
    parse = argparse.ArgumentParser()
    parse.add_argument("--input", help="input sam file", required=True)
    parse.add_argument("--output", help="output filter file", required=True)
    parse.add_argument(
        "--max_indels",
        help="drop the alignment with indels length exceed threshold",
        default=300,
        type=int,
    )
    parse.add_argument("--min_base_quality", default=60, type=int)
    parse.add_argument("--line_num",
                       help="sam file line number",
                       default=None,
                       type=int)
    args = parse.parse_args()
    filter_sam(args.input, args.output, args.max_indels, args.min_base_quality,
               args.line_num)
    end_time = time.time()
    print("alignment filter time cost: %.5f" % (end_time - start_time))
