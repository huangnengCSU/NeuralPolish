bases = ["A", "C", "G", "T", "B"]
label_encoding = {
    "B": 0,
    "<PAD>": 0,
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "<BOS>": 5,
    "<EOS>": 6,
}  # <BOS>:0, <EOS>:5
label_encoding_int2char = {
    0: "<PAD>",
    1: "A",
    2: "C",
    3: "G",
    4: "T",
    5: "<BOS>",
    6: "<EOS>",
}


def kmer_encoding(ksize=5):
    encode_dict = {}
    prev_kmers = [""]
    later_kmers = []
    k = ksize

    while k > 0:
        for v in prev_kmers:
            for i in range(5):
                later_kmers.append(v + bases[i])
        k -= 1
        prev_kmers = later_kmers
        later_kmers = []

    # remove some kmer
    for v in prev_kmers:
        idx = v.find("B")
        if idx == -1:
            kmer = v
        else:
            kmer = v[:idx]
        encode_dict[kmer] = encode_dict.get(kmer, [0] * ksize)

    # 五进制
    for kmer in encode_dict.keys():
        for i in range(len(kmer)):
            if kmer[i] == "A":
                encode_dict[kmer][ksize - len(kmer) + i] = 1
            elif kmer[i] == "C":
                encode_dict[kmer][ksize - len(kmer) + i] = 2
            elif kmer[i] == "G":
                encode_dict[kmer][ksize - len(kmer) + i] = 3
            elif kmer[i] == "T":
                encode_dict[kmer][ksize - len(kmer) + i] = 4
    # 转换成十进制
    for kmer, value in encode_dict.items():
        sum = 0
        for i in range(ksize):
            sum += value[i] * pow(5, ksize - 1 - i)
        encode_dict[kmer] = sum
    # 去掉kmer==""
    encode_dict.pop("")

    # 重新按顺序编码
    output_dict = {"B": 0, "D": 1}
    count_num = 2
    for key in sorted(encode_dict.keys()):
        output_dict[key] = count_num
        count_num += 1

    return output_dict, count_num


def kmer_decoding(ksize=5):
    encode_dict, _ = kmer_encoding(ksize)
    decode_dict = {}
    for kmer, value in encode_dict.items():
        decode_dict[value] = kmer
    return decode_dict, len(decode_dict.keys())


kmer_encoder, kmer_encoder_size = kmer_encoding(ksize=5)
kmer_decoder, kmer_decoder_size = kmer_decoding(ksize=5)

# print(kmer_encoder)

# print("kmer_encoder_size:", kmer_encoder_size)
# print("kmer_decoder_size:", kmer_decoder_size)
