"""
 * @Author: huangneng 
 * @Date: 2019-12-04 16:01:51 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2019-12-04 16:01:51 
"""
import os
import shutil
import re
import time
import argparse
import numpy as np
from kmer_coding import kmer_encoder, kmer_decoder, label_encoding

# next_insert_pat = "\+[0-9]+[ACGTNacgtn]+"
# next_delete_pat = "-[0-9]+[ACGTNacgtn]+"
# forward_mismatch_pat = "[ACGTN]"
# reverse_mismatch_pat = "[acgtn]"
# read_start_pat = "\^."
# read_end_pat = "\$"
# forward_map_pat = "\."
# reverse_map_pat = ","
# single_base_missing_pat = "\*"

pattern = "(\+)([0123456789]+)([ACGTNacgtn]+)|(\-)([0123456789]+)([ACGTNacgtn]+)|\^.|\$|\.|,|[ACGTN]|[acgtn]|\*"

MIN_COVERAGE = 20


class base_structure:
    def __init__(self, current_base, insert_base):
        self.current_base = current_base
        self.insert_base = insert_base


class PileupsFeature:
    def __init__(self):
        self.interval_start = None
        self.contig_name = None
        self.query_features = {}
        self.contig_features = []
        self.max_step = 0

    def clear(self):
        self.interval_start = None
        self.contig_name = None
        self.query_features = {}
        self.contig_features = []
        self.max_step = 0

    def to_array(self, with_contig=True, mode="feature"):
        if mode == "feature":
            array = []
            if with_contig:
                numeric_list = []
                for base in self.contig_features:
                    try:
                        # base stru 20191225
                        numeric_list.append(kmer_encoder[base.current_base])
                        # numeric_list.append(kmer_encoder[base])
                    except:
                        # contig中有碱基N的时候，用填充字符来替换
                        numeric_list.append(kmer_encoder["B"])
                    # base stru 20191225
                    # 插入部分
                    numeric_list.append(kmer_encoder[base.insert_base])

                array.append(numeric_list)

            for key in self.query_features.keys():
                numeric_list = []
                for base in self.query_features[key]:
                    # base stru 20191225
                    numeric_list.append(kmer_encoder[base.current_base])
                    numeric_list.append(kmer_encoder[base.insert_base])
                    # numeric_list.append(kmer_encoder[base])
                array.append(numeric_list)
            array = np.array(array)  # depth,length
            assert array.ndim == 2
        else:
            """
            生成label矩阵时用到的
            """
            array = []

            for key in self.query_features.keys():
                numeric_list = []
                for base in self.query_features[key]:
                    numeric_list.append(kmer_encoder[base.current_base])
                    if base.insert_base != "B":
                        numeric_list.append(kmer_encoder[base.insert_base])
                array.append(numeric_list)
            array = np.array(array)
        return array

    def to_array_copy_contig(self):
        array = []
        numeric_list = []
        for base in self.contig_features:
            try:
                # base stru 20191225
                numeric_list.append(kmer_encoder[base.current_base])
                # numeric_list.append(kmer_encoder[base])
            except:
                # contig中有碱基N的时候，用填充字符来替换
                numeric_list.append(kmer_encoder["B"])
            # base stru 20191225
            # 插入部分
            numeric_list.append(kmer_encoder[base.insert_base])
        array.append(numeric_list)
        for _ in range(MIN_COVERAGE):
            array.append(numeric_list)
        array = np.array(array)  # depth,length
        assert array.ndim == 2
        return array

    def to_file(self, array, output_dir):
        np.save(
            output_dir + "/" + self.contig_name + ":" +
            str(self.interval_start) + "_" +
            str(self.interval_start + self.max_step - 1) + ".feature",
            array,
        )


def parse_mPileup(pileupFile,
                  output_dir,
                  win_size=64,
                  max_insertion_size=5,
                  mode="feature"):
    plp_feat = PileupsFeature()
    set_interval = True
    # interval_start = 0
    read_name_list = []
    pre_read_name_lst = []
    with open(pileupFile) as Fplp:
        coverage_lst = []
        for alignment in Fplp:
            alignment_attrs = alignment.strip().split("\t")
            # [
            #     ref_name,
            #     ref_pos,
            #     ref_base,
            #     mapping_depth,
            #     align_seq,
            #     mapping_quality_seq,
            #     read_name_seq,
            # ] = alignment_attrs
            ref_name = alignment_attrs[0]
            ref_pos = int(alignment_attrs[1])
            ref_base = str.upper(alignment_attrs[2])
            mapping_depth = int(alignment_attrs[3])
            align_seq = alignment_attrs[4]
            mapping_quality_seq = alignment_attrs[5]
            read_name_seq = alignment_attrs[6]
            if set_interval:
                plp_feat.interval_start = ref_pos
                set_interval = False
                plp_feat.contig_name = ref_name
                pre_pos = ref_pos - 1
            if ref_pos - pre_pos != 1:
                """
                ref_pos不连续了，没有reads比对到该区域上，应断开
                000134F	5553	A	0	*	*	*
                000134F	5554	C	0	*	*	*
                000134F	5555	A	0	*	*	*
                000134F	5556	G	0	*	*	*
                000134F	5557	G	1	,$	4	ERR2173373.265096
                000134F	12191	G	1	^",	;	ERR2173373.17992
                000134F	12192	A	1	,	;	ERR2173373.17992
                000134F	12193	T	1	,	8	ERR2173373.17992
                000134F	12194	A	1	,	9	ERR2173373.17992
                """
                print(f"{pre_pos},{ref_pos} mpileup disconnect")
                plp_feat.clear()
                set_interval = True
                coverage_lst = []
                continue
            if ref_base == "N" or ref_base == "n":
                if mode =="feature":
                    ref_base="B"
                else:
                # """
                # 参考序列上为N，断开contig
                # """
                    print(f"N in the contig position: {ref_pos}, disconnect")
                    plp_feat.clear()
                    set_interval = True
                    coverage_lst = []
                    continue

            pre_pos = ref_pos
            coverage_lst.append(mapping_depth)
            base_lst = parse_align_seq(align_seq, ref_base, max_insertion_size)
            read_name_list = read_name_seq.split(",")

            try:
                assert len(base_lst) == mapping_depth
            except:
                """
                # len(base_lst) is one, mapping_depth is zero  base_lst = ['D']

                chromosome	3417501	A	6	,..,,,	9?;.;9	7a608675-7206-4577-9013-0d4be32ace11,2bc05289-bd97-48b0-8c14-de5127af6b21,5f322630-c5ec-4c0c-b62b-beac8d8425d8,3a0a98c4-bcdb-48f5-9d05-2b79554423f9,843a76ce-5339-41bb-94ea-01a27446a026,e54bf263-b71e-4f16-9f8f-5b6246af938f
                chromosome	3417502	A	4	...,	>;:7	2bc05289-bd97-48b0-8c14-de5127af6b21,5f322630-c5ec-4c0c-b62b-beac8d8425d8,0e893b97-896f-4cbb-924d-1e2eba0e230c,e54bf263-b71e-4f16-9f8f-5b6246af938f
                chromosome	3417503	G	0	*	*	*
                chromosome	3417504	C	4	..t,	=.51	c45f0e25-52f2-44fc-947f-eaeab6d252b1,dfad8e12-7f8a-41ef-aaf8-8351cc8617f8,fe58a8f6-d3f0-4b26-ab43-ebf3a7f93fd8,1b0924d3-f0c8-4b63-98c4-14efd9264a40
                chromosome	3417505	T	12	.,*,.,,,,,,,	>;22/=A4H=4;	c45f0e25-52f2-44fc-947f-eaeab6d252b1,7a608675-7206-4577-9013-0d4be32ace11,2bc05289-bd97-48b0-8c14-de5127af6b21,c6e8db5e-7739-4397-9618-c94b847ad1cc,dfad8e12-7f8a-41ef-aaf8-8351cc8617f8,3a0a98c4-bcdb-48f5-9d05-2b79554423f9,68aea4c6-9608-4dc6-b0d1-318edeab7e18,843a76ce-5339-41bb-94ea-01a27446a026,fe58a8f6-d3f0-4b26-ab43-ebf3a7f93fd8,1b0924d3-f0c8-4b63-98c4-14efd9264a40,ba8f2ad1-3fda-44a1-80b3-85d8151ab86f,e54bf263-b71e-4f16-9f8f-5b6246af938f

                # 在3417503位置，没有一条reads比对上，该位置左侧所有比对结果在3417503位置处结束，该位置右侧所有比对结果在3417503位置处重新开始
                # 两边都有比对，中间某个位置的覆盖度为0，是因为samtools mpileup的base quality默认为13，该处的碱基被过滤掉了，如果参数设置为0，则会出现。在这里，我们用参考序列的碱基来作为该位置比对上的碱基，read name因为为空，所以沿用前面的read_name_list
                # -----x------
                # aaaaa*bbbbbb
                # aaaaa*bbbbbb
                # a...a*bbbbbb
                """
                # base stru 20191225
                # if mapping_depth == 0 and len(base_lst) != mapping_depth:
                #     base_lst = [base_structure(ref_base, "B")]
                #     mapping_depth = 1
                # if mapping_depth == 0:
                #     plp_feat.clear()
                #     set_interval = True
                #     coverage_lst = []
                #     continue
                # if base_lst == ["D"] and mapping_depth == 0:
                #     mapping_depth += 1

                # 沿用pre_read_name_lst,用contig上的碱基进行代替, 如果pre_read_name_lst为空，则跳过
                if mapping_depth == 0:
                    # print(
                    #     "mapping_depth is zero, chromosome %s ref pos %d pre_read_name_lst size is %d"
                    #     % (ref_name, ref_pos, len(pre_read_name_lst)))
                    if len(pre_read_name_lst) != 0:
                        base_lst = [base_structure(ref_base, "B")
                                    ] * len(pre_read_name_lst)
                        mapping_depth = len(pre_read_name_lst)
                        read_name_list = pre_read_name_lst
                    else:
                        ## 如果没有pre_read_name_list，就用contig数据，read命名为“extension”
                        base_lst = [base_structure(ref_base, "B")]
                        mapping_depth = 1
                        read_name_list = ["extension"]
                        # plp_feat.clear()
                        # set_interval = True
                        # coverage_lst = []
                        # continue
                """
                ##特殊case1:

                bctg00000005	1745957	C	1	.	~
                bctg00000005	1745958	A	1	.	~
                bctg00000005	1745959	A	1	.	~
                bctg00000005	1745960	A	1	M	~
                bctg00000005	1745961	A	1	.	~
                bctg00000005	1745962	C	1	.	~
                bctg00000005	1745963	A	1	.	~
                # # base M is not in [ACGTNacgtn]
                # # parse_align_seq output:  base_lst = [], len(base_lst)=0, depth = 1
                """
                """
                ##特殊case2:

                000020F 3021720 T       1       .       ~       NC_003074.8
                000020F 3021721 A       1       .       ~       NC_003074.8
                000020F 3021722 T       1       .       ~       NC_003074.8
                000020F 3021723 A       1       .       ~       NC_003074.8
                000020F 3021724 C       1       .       ~       NC_003074.8
                000020F 3021725 T       1       .+30WAAAAAAAAAAAAAAAAAAAAAAAAAAAAA      ~       NC_003074.8
                000020F 3021726 A       1       .       ~       NC_003074.8
                000020F 3021727 A       1       .       ~       NC_003074.8
                000020F 3021728 A       1       .       ~       NC_003074.8
                000020F 3021729 A       1       .       ~       NC_003074.8
                000020F 3021730 A       1       .       ~       NC_003074.8
                # # base W is not in [ACGTNacgtn]
                # # parse_align_seq output:  base_lst = [WAAAAAAAAAAAAAAAAAAAAAAAAAAAAA], len(base_lst)=30, depth = 1
                """
                # base stru 20191225
                if mapping_depth > 0 and len(base_lst) != mapping_depth:
                    base_lst = [base_structure(ref_base, "B")] * mapping_depth
                # if mapping_depth == 1 and base_lst == []:
                #     base_lst.append(ref_base)

                assert len(base_lst) == mapping_depth
            assert len(base_lst) == len(read_name_list)
            read_appear_count = {}
            for i in range(len(read_name_list)):
                qname = read_name_list[i]
                qbase = base_lst[i]
                read_appear_count[qname] = read_appear_count.get(qname, 0) + 1
                if read_appear_count[qname] > 1:
                    continue
                plp_feat.query_features[qname] = plp_feat.query_features.get(
                    qname, [])
                current_feature_len = len(plp_feat.query_features[qname])
                pad_size = plp_feat.max_step - current_feature_len
                if pad_size > 0:
                    # base stru 20191225
                    for _ in range(pad_size):
                        plp_feat.query_features[qname].append(
                            base_structure("B", "B"))
                    # plp_feat.query_features[qname].extend(["B"] * pad_size)
                # if current_feature_len < plp_feat.max_step:
                #     for _ in range(plp_feat.max_step - current_feature_len):
                #         plp_feat.query_features[qname].append("B")
                plp_feat.query_features[qname].append(qbase)

            # 记录 pre_read_name_lst
            if len(read_name_list) > 0:
                pre_read_name_lst = read_name_list

            plp_feat.max_step += 1

            # contig序列
            # base stru 20191225
            plp_feat.contig_features.append(base_structure(ref_base, "B"))
            # plp_feat.contig_features.append(ref_base)

            if ref_pos - plp_feat.interval_start + 1 == win_size:
                for key in plp_feat.query_features.keys():
                    pad_size = plp_feat.max_step - \
                        len(plp_feat.query_features[key])
                    # base stru 20191225
                    for _ in range(pad_size):
                        plp_feat.query_features[key].append(
                            base_structure("B", "B"))
                    # plp_feat.query_features[key].extend(["B"] * (pad_size))
                for key in plp_feat.query_features.keys():
                    assert len(
                        plp_feat.query_features[key]) == plp_feat.max_step
                if mode == "feature":
                    feat_array = plp_feat.to_array(with_contig=True)
                    # if np.mean(coverage_lst) > 0.75 * MIN_COVERAGE:
                    #     feat_array = plp_feat.to_array()
                    # else:
                    #     ### contig比对上的read非常少，低于阈值，则用racon的feature进行填充
                    #     feat_array = plp_feat.to_array_copy_contig()
                elif mode == "label":
                    feat_array = plp_feat.to_array(with_contig=False,
                                                   mode=mode)  # [depth,length]
                    if feat_array.shape[0] > 1:
                        # reference比对到contig上，如果某些contig的覆盖度超过1，则丢弃
                        plp_feat.clear()
                        set_interval = True
                        coverage_lst = []
                        continue
                        # ### reference比对到contig上，有些位置的覆盖度超过1
                        # feat_array = np.expand_dims(feat_array[0], axis=0)
                else:
                    os._exit()
                plp_feat.to_file(feat_array, output_dir)
                plp_feat.clear()
                set_interval = True
                coverage_lst = []
    for key in plp_feat.query_features.keys():
        pad_size = plp_feat.max_step - len(plp_feat.query_features[key])
        for _ in range(pad_size):
            plp_feat.query_features[key].append(base_structure("B", "B"))
        # plp_feat.query_features[key].extend(["B"] * (pad_size))
    for key in plp_feat.query_features.keys():
        assert len(plp_feat.query_features[key]) == plp_feat.max_step
    # print(plp_feat.query_features)
    return


def parse_align_seq(align_seq, ref_base, max_insertion_size):
    """
    解析比对结果中某个位置的所有碱基
    """
    base_list = []
    it = re.finditer(pattern, align_seq)
    for match in it:
        base = match.group()
        if base.startswith("^") or base == "$":
            continue
        elif base == "." or base == ",":
            # base stru 20191225
            base_stru = base_structure(ref_base, "B")
            base_list.append(base_stru)
            ##
            # base_list.append(ref_base)
        elif base == "*":
            # base stru 20191225
            base_stru = base_structure("D", "B")
            base_list.append(base_stru)
            ##
            # base_list.append("D")
        elif base.startswith("+"):
            mark = match.group(1)
            indel_len = match.group(2)
            indel_base = match.group(3)
            if len(indel_base) == indel_len:
                # base_list[-1] += base
                # base stru 20191225
                # 插入错误是后出现的，需要将insert_base加入到前一个base_stru里面
                base_list[-1].insert_base = str.upper(indel_base).replace(
                    "N", "")[:max_insertion_size]
                # 如果出现+3GNC这种情况，N需要被去掉
                # base_list[-1] += str.upper(indel_base).replace("N", "")[
                #     : max_insertion_size - 1
                # ]
            else:
                # 插入序列和mismatch相连的case +3GATct
                # base_list[-1] += "".join(
                #     [mark, indel_len, indel_base[: int(indel_len)]]
                # )
                # insertion的碱基加入到当前位置上,不需要符号和长度
                if base_list[-1] != "D" and base_list[-1] != "B":
                    # base stru 20191225
                    base_list[-1].insert_base = str.upper(
                        indel_base[:int(indel_len)]).replace(
                            "N", "")[:max_insertion_size]
                    # 如果出现GNC这种情况，N需要被去掉
                    # base_list[-1] += str.upper(indel_base[: int(indel_len)]).replace(
                    #     "N", ""
                    # )[: max_insertion_size - 1]
                else:
                    # base stru 20191225
                    base_list[-1].insert_base = str.upper(
                        indel_base[:int(indel_len)]).replace(
                            "N", "")[:max_insertion_size]
                    # "DTTTT" -> "TTTT" and "BTTTT" -> "TTTT"
                    # 如果出现GNC这种情况，N需要被去掉
                    # base_list[-1] = str.upper(indel_base[: int(indel_len)]).replace(
                    #     "N", ""
                    # )[: max_insertion_size - 1]
                for char in indel_base[int(indel_len):]:
                    # 超出indel长度之外的碱基，是紧跟在indel后面的mismatch，与indel无关
                    if char == "N" or char == "n":
                        # base stru 20191225
                        base_stru = base_structure(ref_base, "B")
                        base_list.append(base_stru)
                        # 出现+3GTTNC这种情况，mismatch N用ref_base替换
                        # base_list.append(ref_base)
                    else:
                        # base stru 20191225
                        base_stru = base_structure(str.upper(char), "B")
                        base_list.append(base_stru)
                        # base_list.append(str.upper(char))
        elif base.startswith("-"):
            mark = match.group(4)
            indel_len = match.group(5)
            indel_base = match.group(6)
            if len(indel_base) == indel_len:
                # base_list[-1] += base
                pass
            else:
                # 删除序列和mismatch相连的case -3GATct
                # base_list[-1] += "".join(
                #     [mark, indel_len, indel_base[: int(indel_len)]]
                # )
                # deletion的碱基不需要加入到当前位置上,也不需要符号和长度
                for char in indel_base[int(indel_len):]:
                    # 超出indel长度之外的碱基，是紧跟在indel后面的mismatch，与indel无关
                    if char == "N" or char == "n":
                        # base stru 20191225
                        base_stru = base_structure(ref_base, "B")
                        base_list.append(base_stru)
                        # 出现-3GTTNC这种情况，mismatch N用ref_base替换
                        # base_list.append(ref_base)
                    else:
                        # base stru 20191225
                        base_stru = base_structure(str.upper(char), "B")
                        base_list.append(base_stru)
                        # base_list.append(str.upper(char))
        else:
            # mismatch
            if base == "N" or base == "n":
                # reference has base N
                # bctg00000005	2502779	C	1	N	~
                # bctg00000005	2502780	A	1	N	~
                # bctg00000005	2502781	C	1	N	~
                # bctg00000005	2502782	A	1	.	~
                # bctg00000005	2502783	C	1	.	~
                # bctg00000005	2502784	A	1	.	~
                # bctg00000005	2502785	C	1	.-12ATACACAATACA	~
                # bctg00000005	2502786	A	1	*	~
                # bctg00000005	2502787	T	1	*	~
                # bctg00000005	2502788	A	1	*	~
                # bctg00000005	2502789	C	1	*	~
                ###

                # base stru 20191225
                base_stru = base_structure(ref_base, "B")
                base_list.append(base_stru)
                # base_list.append(ref_base)
            else:
                # base stru 20191225
                base_stru = base_structure(str.upper(base), "B")
                base_list.append(base_stru)
                # base_list.append(str.upper(base))
    return base_list


if __name__ == "__main__":
    start_time = time.time()
    parse = argparse.ArgumentParser()
    parse.add_argument("-i", help="input the mpileup file", required=True)
    parse.add_argument("-o", help="output dir", required=True)
    parse.add_argument("-w", help="window size", type=int, default=64)
    parse.add_argument("-m", help="max insertion size", type=int, default=5)

    opt = parse.parse_args()
    if not os.path.exists(opt.o):
        os.mkdir(opt.o)
        os.mkdir(opt.o + "/features")
        os.mkdir(opt.o + "/laebls")
    else:
        rm_confirm = "no"
        rm_confirm = input(
            "do you want to delete the directory %s \n please input `yes` or `no`:"
            % opt.o)
        if rm_confirm == "yes":
            shutil.rmtree(opt.o)
            print("remove directory %s completed." % opt.o)
            os.mkdir(opt.o)
            os.mkdir(opt.o + "/features")
            os.mkdir(opt.o + "/laebls")
    parse_mPileup(opt.i, opt.o + "/features", opt.w, opt.m)
    print("total time cost: %.5f min" % ((time.time() - start_time) / 60))

# st = time.time()
# parse_mPileup(
#     "/mnt/c/Users/hn/Desktop/mpileup_bctg00000000_1_10000.txt",
#     "/mnt/h/github/MDPileups/output_data",
# )
# print(time.time() - st)
