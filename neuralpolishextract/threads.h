//
// Created by csuhu on 2020/7/17.
//

#ifndef EXTRACT_THREADS_H
#define EXTRACT_THREADS_H

#endif //EXTRACT_THREADS_H

#include "load_pileup.h"
#include <thread>
#include <mutex>

template<typename W>
void MultiThreadRun(uint32_t thread_size, W work_func);

void
generate_features(vector<string> &filelist, uint32_t window_length, map<string, int> &encoder, uint32_t thread_size,
                  uint32_t max_coverage, string &outdir);
