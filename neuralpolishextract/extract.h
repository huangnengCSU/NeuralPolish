//
// Created by lab on 2020/7/17.
//

#ifndef NEURALPOLISHEXTRACT_EXTRACT_H
#define NEURALPOLISHEXTRACT_EXTRACT_H

#endif //NEURALPOLISHEXTRACT_EXTRACT_H


#include "threads.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>


bool get_all_files(const string& dir_in, vector<string>& files);
std::map<std::string, int> load_encoder(const char *encoder_file);
std::map<int, std::string> load_decoder(const char *decoder_file);



