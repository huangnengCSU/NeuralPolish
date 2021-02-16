//
// Created by csuhu on 2020/7/17.
//

#include "threads.h"

template<typename W>
void MultiThreadRun(uint32_t thread_size, W work_func) {
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(work_func, i));
    }

    for (uint32_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}

void
generate_features(vector<string> &filelist, uint32_t window_length, map<string, int> &encoder, uint32_t thread_size,
                  uint32_t max_coverage, string &outdir) {
    std::mutex mutex_gen;
    std::mutex mutex_comb;
    vector<string>::iterator fileiter = filelist.begin();
    auto generate_func = [&mutex_gen, &filelist, &fileiter]() {
        std::lock_guard<std::mutex> lock(mutex_gen);

        int get_count = 0;
        string filename;
        if (fileiter != filelist.end()) {
            filename = *fileiter;
            fileiter++;
            get_count++;
        }
        return std::pair<int, string>(get_count, filename);
    };

    auto combine_func = [&mutex_comb](vector<std::pair<string, Matrix>> &mat_pool, string outdir) {
        std::lock_guard<std::mutex> lock(mutex_comb);
        for (auto item : mat_pool) {
            string outfilename = outdir + "/" + item.first + ".feature";
            item.second.to_file(outfilename);
        }
    };

    auto work_func = [generate_func, combine_func, &filelist, &encoder, &window_length, &max_coverage, &outdir](
            size_t) {
        std::vector<std::pair<std::string, Matrix>> mat_pool;
        while (true) {
            std::pair<int, string> product = generate_func();
            pair<vector<string>, vector<Matrix>> result;
            if (product.first > 0) {
                result = read_pileup(product.second, window_length, max_coverage, encoder);
                for (int i = 0; i < result.first.size(); i++) {
                    mat_pool.push_back(std::pair<std::string, Matrix>(result.first[i], result.second[i]));
                }
                if (mat_pool.size() > 5000) {
                    combine_func(mat_pool, outdir);
                    mat_pool.clear();
                }
            } else {
                break;
            }
        }
        if (mat_pool.size() > 0) {
            combine_func(mat_pool, outdir);
            mat_pool.clear();
        }
    };
    MultiThreadRun(thread_size, work_func);
}