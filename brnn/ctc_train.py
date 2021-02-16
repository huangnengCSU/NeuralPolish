'''
 * @Author: huangneng 
 * @Date: 2020-01-19 14:16:05 
 * @Last Modified by:   huangneng 
 * @Last Modified time: 2020-01-19 14:16:05 
 '''

import argparse
import yaml
from rnn_ctc import train, Model, AttrDict
from optim import Optimizer

import torch
from torch.optim import Adam

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', '-c', type=str, required=True)
    parser.add_argument('--train_data', '-t', required=True)
    parser.add_argument('--eval_data', '-e', required=True)
    parser.add_argument('--no_cuda', action="store_true")
    opt = parser.parse_args()

    configfile = open(opt.config)
    config = AttrDict(yaml.load(configfile))
    device = torch.device('cuda' if not opt.no_cuda else 'cpu')
    config.device = device
    config.train_data = opt.train_data
    config.eval_data = opt.eval_data
    model = Model(config.Depth, config.Vocab_size, config.model.Embed_size,
                  config.model.Rnn1_hid, config.model.Rnn1_nlayer,
                  config.model.Rnn1_dropout, config.model.Hid_dim,
                  config.model.Rnn2_hid, config.model.Rnn2_nlayer,
                  config.model.Rnn2_dropout,
                  config.model.Out_dim).to(config.device)
    # optim = Adam(model.parameters(), config.optim.lr)
    optim = Optimizer(model.parameters(), config.optim)
    train(model, optim, config)
