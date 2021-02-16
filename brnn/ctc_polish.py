'''
 * @Author: huangneng 
 * @Date: 2020-01-23 16:35:16 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2020-01-23 16:35:16 
 '''

import argparse
import yaml
from rnn_ctc import Model, polish, AttrDict
import torch
import time

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("-config", type=str, required=True)
    parser.add_argument("-model_path", required=True)
    parser.add_argument('-data_path', required=True)
    parser.add_argument('-prepolished_path',required=True)
    parser.add_argument('-output', required=True)
    parser.add_argument('--no_cuda', action="store_true")
    opt = parser.parse_args()

    configfile = open(opt.config)
    config = AttrDict(yaml.load(configfile))
    device = torch.device('cuda' if not opt.no_cuda else 'cpu')
    config.device = device
    model = Model(config.Depth, config.Vocab_size, config.model.Embed_size,
                  config.model.Rnn1_hid, config.model.Rnn1_nlayer,
                  config.model.Rnn1_dropout, config.model.Hid_dim,
                  config.model.Rnn2_hid, config.model.Rnn2_nlayer,
                  config.model.Rnn2_dropout,
                  config.model.Out_dim).to(config.device)
    polish(model, opt.model_path, opt.data_path, opt.prepolished_path, opt.output, config)
    end_time = time.time()
    print("total inference time cost: %.5f" % (end_time - start_time))
