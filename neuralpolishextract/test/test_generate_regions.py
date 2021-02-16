def get_mPileup_regions(ctg_len, win_size):
    mpileup_region_list = []
    contig = "contig_1"
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

if __name__ == '__main__':
    regs = get_mPileup_regions(5000,10000)
    for reg in regs:
        print(reg)
    print()
    regs = get_mPileup_regions(9999,10000)
    for reg in regs:
        print(reg)
    print()
    regs = get_mPileup_regions(10000,10000)
    for reg in regs:
        print(reg)
    print()
    regs = get_mPileup_regions(10001,10000)
    for reg in regs:
        print(reg)
    print()
    regs = get_mPileup_regions(10050,10000)
    for reg in regs:
        print(reg)
    print()
    regs = get_mPileup_regions(50000,10000)
    for reg in regs:
        print(reg)
    print()
