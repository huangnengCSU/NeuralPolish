//
// Created by lab on 2020/7/17.
//

#ifndef EXTRACT_LOAD_PILEUP_H
#define EXTRACT_LOAwhileD_PILEUP_H

#endif //EXTRACT_LOAD_PILEUP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <assert.h>
#include <regex>
#include <string>
#include <algorithm>

using namespace std;

struct plp_elem {
    plp_elem(string cur, string ins) {
        this->current_nucle = cur;
        this->insertion_nucles = ins;
    }

    string current_nucle;
    string insertion_nucles;
};

struct matrix_ele {
    matrix_ele(int cur, int ins) {
        this->current_nucle = cur;
        this->insertion_nucles = ins;
    }

    int current_nucle;
    int insertion_nucles;
};

class Matrix {
public:
    int row_num = 0, col_num = 0;
    std::vector<std::vector<matrix_ele>> m;

public:
    Matrix() {};

    Matrix(int max_coverage);

    ~Matrix() {};

    void Insert(matrix_ele init_elem);

    void to_file(string filename);
};


void SplitString(string &s, vector<string> &v, string &c);

pair<vector<string>, vector<Matrix>>
read_pileup(string filename, uint32_t window_length, uint32_t max_coverage, map<string, int> &encoder);

void
parse_mapping_bases(string &mapping_bases, regex &r, vector<plp_elem> &parsed_vec, string &ref_base, int max_coverage);

void fill_in_matrix(vector<plp_elem> &align_base_vector, vector<string> &mapping_readnames, Matrix &mat,
                    map<string, int> &encoder, map<string, int> &read_name_idx, int &row_count);
