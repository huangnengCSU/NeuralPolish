"""
 * @Author: huangneng 
 * @Date: 2019-12-09 00:37:04 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2019-12-09 00:37:04 
 """
__author__ = "hn"

import argparse
import os
import random
import shutil

parse = argparse.ArgumentParser()
parse.add_argument("-train", required=True)
parse.add_argument("-eval", required=True)
parse.add_argument("-p", help="validate ratio", default=0.1, type=float)
parse.add_argument("-choose_count", help="count of validate", default=None, type=int)
parse.add_argument("-copy", action="store_true", default=False)
args = parse.parse_args()

if not os.path.exists(args.eval):
    os.mkdir(args.eval)
    os.mkdir(os.path.join(args.eval, "features"))
    os.mkdir(os.path.join(args.eval, "labels"))
else:
    rm_confirm = "no"
    rm_confirm = input(
        "do you want to delete the directory %s \n please input `yes` or `no`:"
        % args.eval
    )
    if rm_confirm == "yes":
        shutil.rmtree(args.eval)
        print("remove directory %s completed." % args.eval)
        os.mkdir(args.eval)
        os.mkdir(os.path.join(args.eval, "features"))
        os.mkdir(os.path.join(args.eval, "labels"))


signal_name_lst = os.listdir(args.train + "/features/")
label_name_lst = os.listdir(args.train + "/labels/")
intersection = list(set(signal_name_lst) & set(label_name_lst))
size = len(intersection)
if args.choose_count is None:
    choose_size = int(size * args.p)
else:
    choose_size = args.choose_count

choose_names = random.sample(intersection, choose_size)

for v in choose_names:
    if not args.copy:
        shutil.move(args.train + "/features/" + v, args.eval + "/features/")
    else:
        shutil.copy(args.train + "/features/" + v, args.eval + "/features/")

for v in choose_names:
    if not args.copy:
        shutil.move(args.train + "/labels/" + v, args.eval + "/labels/")
    else:
        shutil.copy(args.train + "/labels/" + v, args.eval + "/labels/")

assert len(os.listdir(args.eval + "/features")) == len(
    os.listdir(args.eval + "/labels")
)
