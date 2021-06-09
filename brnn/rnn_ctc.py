"""
 * @Author: huangneng
 * @Date: 2020-01-17 16:12:47
 * @Last Modified by:   huang neng
 * @Last Modified time: 2021-06-09 10:26:49
 """

import torch
import torch.nn as nn
import torch.utils.data as Data
import torch.nn.functional as F
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence
from tensorboardX import SummaryWriter
from multiprocessing import Pool

import os
import math
import random
import numpy as np

from ctcdecoder import GreedyDecoder
from kmer_coding import kmer_decoder, kmer_decoder_size

"""  模型结构
input region (Batch,length,20)
|
embedding layer(embed_size)
|
(Batch,length,20,embed_size)
|
(Batch*20,length,embed_size)
|
rnn layer (inside each read)
|
(Batch*20,length,hid1)
|
(Batch,length,20*hid1)
|
Linear (20*hid1,hid2)
|
(Batch,length,hid2)
|
rnn layer (between reads)
|
(Batch,length,output_dim)
|
ctc layer
|
output
"""


class TrainDataset(Data.Dataset):
    def __init__(self, data_path):
        super(TrainDataset, self).__init__()
        self.data_path = data_path
        self.features = os.listdir(data_path + "/features")
        self.labels = os.listdir(data_path + "/labels")
        self.intersection = list(set(self.features) & set(self.labels))

    def __len__(self):
        return len(self.intersection)

    def __getitem__(self, i):
        feature_name = self.intersection[i]
        feature = np.load(self.data_path + "/features/" +
                          feature_name)  # [depth,length]
        label = np.squeeze(np.load(self.data_path + "/labels/" + feature_name),
                           axis=0)  # (length,)
        # feature_depth = feature.shape[0]
        # feature_length = feature.shape[1]
        return feature, label, feature_name


class EvalDataset(Data.Dataset):
    def __init__(self, data_path):
        super(EvalDataset, self).__init__()
        self.data_path = data_path
        self.features = os.listdir(data_path + "/features")
        self.labels = os.listdir(data_path + "/labels")
        self.intersection = list(set(self.features) & set(self.labels))

    def __len__(self):
        return len(self.intersection)

    def __getitem__(self, i):
        feature_name = self.intersection[i]
        feature = np.load(self.data_path + "/features/" +
                          feature_name)  # [depth,length]
        label = np.squeeze(np.load(self.data_path + "/labels/" + feature_name),
                           axis=0)  # (length,)
        # feature_depth = feature.shape[0]
        # feature_length = feature.shape[1]
        return feature, label, feature_name


class Collate:
    def __init__(self, sample_depth, random_seed):
        self.sample_depth = sample_depth
        self.seed = random_seed

    def _collate(self, batch_data):
        batchsize = len(batch_data)
        feature_length_list, label_length_list = [], []
        for data_item in batch_data:
            feature_length_list.append(data_item[0].shape[1])
            label_length_list.append(data_item[1].shape[0])
        max_feature_length = max(feature_length_list)
        max_label_length = max(label_length_list)

        feature_lst, label_lst, feature_name_lst = [], [], []
        feature_tensor_len_lst, label_tensor_len_lst = [], []
        for data_item in batch_data:
            feature = torch.LongTensor(data_item[0])  # [depth,length]
            label = torch.LongTensor(data_item[1])  # [length,]
            feature_name = data_item[2]
            feature_name_lst.append(feature_name)

            depth = feature.shape[0]
            feature_length = feature.shape[1]
            label_length = label.shape[0]

            if depth < self.sample_depth:
                expand_ratio = math.ceil(self.sample_depth * 2 / depth)
                expand_feature = torch.cat([feature] * expand_ratio, dim=0)
                random.seed(self.seed)
                sample_index = random.sample(list(range(depth * expand_ratio)),
                                             self.sample_depth)
                sample_feature = expand_feature[sample_index]
            else:
                random.seed(self.seed)
                sample_index = random.sample(list(range(depth)),
                                             self.sample_depth)
                sample_feature = feature[sample_index]
            sample_feature = F.pad(
                sample_feature, [0, max_feature_length - feature_length, 0, 0],
                "constant",
                value=0)
            label = F.pad(label, [0, max_label_length - label_length],
                          "constant",
                          value=0)
            feature_lst.append(sample_feature)
            feature_tensor_len_lst.append(feature_length)
            label_lst.append(label)
            label_tensor_len_lst.append(label_length)
        feature_tensor = torch.cat(feature_lst,
                                   0).view(batchsize, self.sample_depth,
                                           -1).contiguous().transpose(1, 2)
        label_tensor = torch.cat(label_lst, 0).view(batchsize, -1)
        feature_len_tensor = torch.LongTensor(feature_tensor_len_lst)
        label_len_tensor = torch.LongTensor(label_tensor_len_lst)

        # # 按照序列长度从大到小排序
        # feature_len_tensor, perm_idx = feature_len_tensor.sort(0,
        #                                                        descending=True)
        # feature_tensor = feature_tensor[perm_idx]
        # label_tensor = label_tensor[perm_idx]
        # label_len_tensor = label_len_tensor[perm_idx]

        return (feature_tensor, label_tensor, feature_len_tensor,
                label_len_tensor, feature_name_lst)

    def __call__(self, batch_data):
        return self._collate(batch_data)


class Feature_Store:
    feature = None
    name = None
    start_pos = None

    def __init__(self, feature, name, pos):
        self.feature = feature
        self.name = name
        self.start_pos = pos


def load_data(fpath):
    data_array = []
    with open(fpath, 'r') as fin:
        for line in fin:
            line_lst = [int(v) for v in line.strip().split(',')]
            data_array.append(line_lst)
    data_array = np.array(data_array)
    if data_array.ndim == 1:
        data_array.reshape((1, -1))
    return data_array


class PolishDataset(Data.Dataset):
    def __init__(self, contig_feature_paths):
        super(PolishDataset, self).__init__()
        self.data_path = os.path.dirname(contig_feature_paths[0])
        data_files = [os.path.basename(v) for v in contig_feature_paths]
        self.features = sorted(
            data_files,
            key=lambda x: int(x.split(":")[1].split("-")[0]),
            reverse=False)
        self.count = len(self.features)
        print(os.path.basename(contig_feature_paths[0]).split(':')[0], '\t', self.count)

    def __len__(self):
        return self.count

    def __getitem__(self, i):
        feature_name = self.features[i]
        feature = load_data(self.data_path + os.path.sep + feature_name)
        return feature, feature_name


class PolishCollate:
    def __init__(self, sample_depth, random_seed):
        self.sample_depth = sample_depth
        self.seed = random_seed

    def _collate(self, batch_data):
        batchsize = len(batch_data)
        feature_length_list = []
        for data_item in batch_data:
            feature_length_list.append(data_item[0].shape[1])
        max_feature_length = max(feature_length_list)

        feature_lst, feature_name_lst = [], []
        feature_tensor_len_lst = []
        for data_item in batch_data:
            feature = torch.LongTensor(data_item[0])  # [depth,length]
            feature_name = data_item[1]
            feature_name_lst.append(feature_name)

            depth = feature.shape[0]
            feature_length = feature.shape[1]

            if depth < self.sample_depth:
                expand_ratio = math.ceil(self.sample_depth * 2 / depth)
                expand_feature = torch.cat([feature] * expand_ratio, dim=0)
                random.seed(self.seed)
                sample_index = random.sample(list(range(depth * expand_ratio)),
                                             self.sample_depth)
                sample_feature = expand_feature[sample_index]
            else:
                random.seed(self.seed)
                sample_index = random.sample(list(range(depth)),
                                             self.sample_depth)
                sample_feature = feature[sample_index]
            sample_feature = F.pad(
                sample_feature, [0, max_feature_length - feature_length, 0, 0],
                "constant",
                value=0)
            feature_lst.append(sample_feature)
            feature_tensor_len_lst.append(feature_length)
        feature_tensor = torch.cat(feature_lst,
                                   0).view(batchsize, self.sample_depth,
                                           -1).contiguous().transpose(1, 2)
        feature_len_tensor = torch.LongTensor(feature_tensor_len_lst)

        return (feature_tensor, feature_len_tensor, feature_name_lst)

    def __call__(self, batch_data):
        return self._collate(batch_data)


class Model(nn.Module):
    def __init__(self, depth, vocab_size, embed_size, rnn1_hid, rnn1_nlayer,
                 rnn1_dropout, hid_dim, rnn2_hid, rnn2_nlayer, rnn2_dropout,
                 out_dim):
        super(Model, self).__init__()
        self.depth = depth
        self.embed_size = embed_size
        self.embedding_layer = nn.Embedding(vocab_size, embed_size)
        self.rnn_layer_inside_read = nn.GRU(embed_size,
                                            rnn1_hid,
                                            rnn1_nlayer,
                                            batch_first=True,
                                            bidirectional=True,
                                            dropout=rnn1_dropout)
        self.linear_between_two_rnns = nn.Linear(depth * rnn1_hid * 2, hid_dim)
        self.rnn_layer_between_reads = nn.GRU(hid_dim,
                                              rnn2_hid,
                                              rnn2_nlayer,
                                              batch_first=True,
                                              bidirectional=True,
                                              dropout=rnn2_dropout)
        self.out_linear = nn.Linear(2 * rnn2_hid, out_dim)

    def forward(self, input_tensor, feature_lengths):
        """
        input_tensor: shape of [batch,length,depth]
        """
        batchsize = input_tensor.size(0)
        length = input_tensor.size(1)

        # [batch,length,depth,embed_size]
        out = self.embedding_layer(input_tensor)
        # [batch,depth,length,embed_size]
        out = out.transpose(1, 2).contiguous()
        # [batch*depth,length,embed_size]
        out = out.view(-1, length, self.embed_size)

        # [batch*depth,length,2*rnn1_hid]
        out, _ = self.rnn_layer_inside_read(out)

        out = out.view(batchsize, self.depth, length, -1)
        # [batch,length,depth,2*rnn1_hid]
        out = out.transpose(1, 2).contiguous()
        # [batch,length,depth*2*rnn1_hid]
        out = out.view(batchsize, length, -1)

        out = self.linear_between_two_rnns(out)  # [batch,length,hid_dim]

        out, _ = self.rnn_layer_between_reads(out)  # [batch,length,2*rnn2_hid]

        out = self.out_linear(out)  # [batch,length,out_dim]

        out = out.log_softmax(dim=2)
        return out


def test_TrainDataset(train_data_path, sample_depth, seed):
    train_dataset = Data.DataLoader(TrainDataset(train_data_path),
                                    batch_size=4,
                                    num_workers=40,
                                    collate_fn=Collate(sample_depth, seed),
                                    shuffle=True)
    for (step, (pileup_features, targets, feature_lengths, target_lengths,
                pileup_file_names)) in enumerate(train_dataset):
        print(pileup_features)
        print(targets)


def test_EvalDataset(eval_data_path, sample_depth, seed):
    eval_dataset = Data.DataLoader(EvalDataset(eval_data_path),
                                   batch_size=4,
                                   num_workers=40,
                                   collate_fn=Collate(sample_depth, seed),
                                   shuffle=True)
    for (step, (pileup_features, targets, feature_lengths, target_lengths,
                pileup_file_names)) in enumerate(eval_dataset):
        print(pileup_features)
        print(targets)


def train(model, optimizer, opt):
    print(model)
    train_dataset = Data.DataLoader(TrainDataset(opt.train_data),
                                    batch_size=opt.train.Batch_size,
                                    num_workers=40,
                                    collate_fn=Collate(opt.Depth, opt.Seed),
                                    shuffle=True)
    batch_steps = len(train_dataset)
    visualizer = SummaryWriter(
        os.path.join(opt.train.exp_name, opt.train.visualizer_dir, "log"))
    train_global_step = 0
    for epoch in range(opt.train.Epochs):
        model.train()
        total_train_loss = 0
        for (step, (pileup_features, targets, feature_lengths, target_lengths,
                    pileup_file_names)) in enumerate(train_dataset):
            train_global_step += step
            pileup_features = pileup_features.to(opt.device)  # [128,128,20]
            targets = targets.to(opt.device)  # [128,65]
            feature_lengths = feature_lengths.to(opt.device)
            target_lengths = target_lengths.to(opt.device)

            optimizer.zero_grad()
            preds = model(pileup_features, feature_lengths)
            preds = preds.transpose(0, 1)  # [length,batch,outdim]
            loss = F.ctc_loss(
                preds,
                targets,
                feature_lengths,
                target_lengths,
                blank=kmer_decoder_size)  # kmer_decoder: 0-1365, blank=1366
            loss.backward()
            total_train_loss += loss.item()
            grad_norm = nn.utils.clip_grad_norm_(model.parameters(),
                                                 opt.train.Clip)
            optimizer.step()
            avg_loss = total_train_loss / (step + 1)
            process = step / batch_steps * 100
            if (step + 1) % opt.train.Show_step == 0:
                print(
                    "-Training-Epoch:%d(%.5f%%), train loss:%.5f, average loss:%.5f"
                    % (epoch, process, loss.item(), avg_loss))
                visualizer.add_scalar("train_loss", loss.item(),
                                      train_global_step)
        eval_loss, eval_cer, eval_step = eval(model, opt)
        print("+Eval-Epoch:%d, average eval loss:%.5f, average eval cer:%.5f" %
              (epoch, eval_loss, eval_cer))

        visualizer.add_scalar("eval_loss", eval_loss, epoch)
        visualizer.add_scalar("eval_cer", eval_cer, epoch)
        save_name = os.path.join(
            opt.train.exp_name,
            "%s.epoch%d.chkpt" % (opt.train.save_model, epoch))
        save_model(model, optimizer, opt, save_name, epoch, train_global_step)
        if epoch >= opt.optim.begin_to_adjust_lr:
            optimizer.decay_lr()
            if optimizer.lr<1e-6:
                print("The learning rate is too low to train.")
                break
            print("Epoch %d update learning rate %.6f"%(epoch,optimizer.lr))
    print("The training process is over.")


def eval(model, opt):
    eval_dataset = Data.DataLoader(EvalDataset(opt.eval_data),
                                   batch_size=opt.train.Batch_size,
                                   num_workers=40,
                                   collate_fn=Collate(opt.Depth, opt.Seed),
                                   shuffle=False)
    decoder = GreedyDecoder(kmer_decoder, blank_index=kmer_decoder_size)
    model.eval()
    total_eval_loss = 0
    total_cer, num_chars = 0, 0
    eval_step_count = 0
    for (step, (pileup_features, targets, feature_lengths, target_lengths,
                pileup_file_names)) in enumerate(eval_dataset):
        eval_step_count += 1
        pileup_features = pileup_features.to(opt.device)  # [128,128,20]
        targets = targets.to(opt.device)  # [128,65]
        feature_lengths = feature_lengths.to(opt.device)
        target_lengths = target_lengths.to(opt.device)

        preds = model(pileup_features, feature_lengths)
        preds = preds.transpose(0, 1)  # [length,batch,outdim]
        loss = F.ctc_loss(preds,
                          targets,
                          feature_lengths,
                          target_lengths,
                          blank=kmer_decoder_size)

        pred_strings, offset = decoder.decode(preds.transpose(0, 1),
                                              feature_lengths)
        target_strings = decoder.convert_to_strings(targets, target_lengths)

        total_eval_loss += loss.item()
        for i in range(len(targets)):
            cer = decoder.cer(pred_strings[i], target_strings[i])
            total_cer += cer
            num_chars += len(target_strings[i])
    average_eval_loss = total_eval_loss / eval_step_count
    average_cer = float(total_cer) / num_chars
    return average_eval_loss, average_cer, eval_step_count


def read_assembly(assembly_path):
    assembly_contigs = {}
    with open(assembly_path, 'r') as fasm:
        for line in fasm:
            header = line[1:].strip()
            contig_name = header.split(' ')[0]
            sequence = fasm.readline().strip()
            assembly_contigs[contig_name] = sequence
    return assembly_contigs


def polish(model, model_path, data_path, racon_assembly_path, output_file, opt):
    assembly_contigs = read_assembly(racon_assembly_path)
    model = load_model(model, model_path, opt).to(opt.device)
    model.eval()
    decoder = GreedyDecoder(kmer_decoder, blank_index=kmer_decoder_size)

    contig_dict = sort_inference_data(data_path)
    consensus_out = open(output_file, "w")
    for contig_name in contig_dict.keys():
        if len(contig_dict[contig_name]) == 1:
            contig_windows_chars = []
            contig_feature_paths = [v[0] for v in contig_dict[contig_name][0]]
            polish_dataset = Data.DataLoader(PolishDataset(contig_feature_paths),
                                             batch_size=opt.train.Batch_size,
                                             num_workers=40,
                                             collate_fn=PolishCollate(
                                                 opt.Depth, opt.Seed),
                                             shuffle=False)
            for (step, (pileup_features, feature_lengths,
                        pileup_filename)) in enumerate(polish_dataset):
                pileup_features = pileup_features.to(opt.device)
                feature_lengths = feature_lengths.to(opt.device)
                preds = model(pileup_features, feature_lengths)
                pred_strings, offset = decoder.decode(preds, feature_lengths)
                contig_windows_chars.extend(pred_strings)
            consensus_out.write(">" + contig_name.replace('$', '_') + "\n")
            for win_chars in contig_windows_chars:
                consensus_out.write(win_chars)
            consensus_out.write("\n")
        else:
            """
            split contigs, fill with assembly
            """
            split_number = len(contig_dict[contig_name])
            consensus_out.write(">" + contig_name.replace('$', '_') + "\n")
            for idx in range(split_number):
                contig_windows_chars = []
                contig_feature_paths = [v[0]
                                        for v in contig_dict[contig_name][idx]]
                polish_dataset = Data.DataLoader(PolishDataset(contig_feature_paths),
                                                 batch_size=opt.train.Batch_size,
                                                 num_workers=40,
                                                 collate_fn=PolishCollate(
                                                     opt.Depth, opt.Seed),
                                                 shuffle=False)
                for (step, (pileup_features, feature_lengths,
                            pileup_filename)) in enumerate(polish_dataset):
                    pileup_features = pileup_features.to(opt.device)
                    feature_lengths = feature_lengths.to(opt.device)
                    preds = model(pileup_features, feature_lengths)
                    pred_strings, offset = decoder.decode(
                        preds, feature_lengths)
                    contig_windows_chars.extend(pred_strings)

                for win_chars in contig_windows_chars:
                    consensus_out.write(win_chars)
                # consensus_out.write("\n")
                if idx < (split_number - 1):
                    # end pos
                    break_left_pos = contig_dict[contig_name][idx][-1][2]
                    # start pos
                    break_right_pos = contig_dict[contig_name][idx + 1][0][1]
                    # end_pos+1:start_pos
                    consensus_out.write(
                        assembly_contigs[contig_name][break_left_pos + 1:break_right_pos])
            consensus_out.write('\n')

    consensus_out.close()


def sort_inference_data(inference_data_path):
    contig_dict = {}
    feature_path = inference_data_path + "/features"
    data_files = os.listdir(feature_path)
    for chk_file in data_files:
        contig_name = chk_file.split(":")[0]
        start_pos, end_pos = chk_file.split(":")[1].split(".")[0].split("-")
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        fpath = feature_path + "/" + chk_file
        contig_dict[contig_name] = contig_dict.get(contig_name, [])
        contig_dict[contig_name].append([fpath, start_pos, end_pos])
        # if contig_dict[contig_name] == 1:
        #     os.mkdir(feature_path + "/" + contig_name)
        # shutil.move(feature_path + "/" + chk_file, feature_path + "/" + contig_name)
    """
    找到break point,把contig断开
    """
    split_contig_dict = {}
    before_break_contig_nums = len(contig_dict.keys())
    breaked_contig_nums = 0
    breaked_points = 0
    for contig_name in contig_dict.keys():
        split_contigs = [[]]
        sorted_contig_pos = sorted(contig_dict[contig_name],
                                   key=lambda v: v[1])
        pre_end_pos = sorted_contig_pos[0][2]
        split_contigs[-1].append(
            sorted_contig_pos[0])  # add path to last list [[000000F:1_64,000000F_65,128],[000000F:1000_1064]]
        for v in sorted_contig_pos[1:]:
            cur_start_pos = v[1]
            cur_end_pos = v[2]
            if cur_start_pos == pre_end_pos + 1:
                split_contigs[-1].append(v)
            else:
                # break point, new contig
                split_contigs.append([])
                split_contigs[-1].append(v)
                breaked_points += 1
            pre_end_pos = cur_end_pos
        split_contig_dict[contig_name] = split_contigs
        # split_num = len(split_contigs)
        # breaked_contig_nums += split_num
        # if split_num == 1:
        #     split_contig_dict[contig_name] = split_contigs[0]
        # else:
        #     for idx in range(split_num):
        #         split_contig_dict[contig_name + "$" +
        #                           str(idx)] = split_contigs[idx]
    # print("before breaking contigs number: ", before_break_contig_nums)
    # print("after breaking contigs number: ", breaked_contig_nums)
    print('bread points number:', breaked_points)
    return split_contig_dict


def save_model(model, optimizer, config, save_name, epoch, step):
    # multi_gpu = True if config.training.num_gpu > 1 else False
    checkpoint = {
        "model": model.state_dict(),
        "optimizer": optimizer.state_dict(),
        "epoch": epoch,
        "step": step,
    }
    torch.save(checkpoint, save_name)


def load_model(model, model_path, opt):
    if opt.device==torch.device('cpu'):
        device = 'cpu'
    else:
        device = 'cuda'
    checkpoint = torch.load(model_path, map_location=device)
    model.load_state_dict(checkpoint["model"])
    return model


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)

    def __getattr__(self, item):
        if item not in self:
            return None
        if type(self[item]) is dict:
            self[item] = AttrDict(self[item])
        return self[item]

    def __setattr__(self, item, value):
        self.__dict__[item] = value


if __name__ == '__main__':
    import sys

    data_path = sys.argv[1]
    sample_depth = int(sys.argv[2])
    seed = int(sys.argv[3])
    test_TrainDataset(data_path, sample_depth, seed)
