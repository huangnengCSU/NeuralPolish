//
// Created by lab on 2020/7/17.
//

#ifndef EXTRACT_TEST_PARSE_MAPPING_BASES_H
#define EXTRACT_TEST_PARSE_MAPPING_BASES_H

#endif //EXTRACT_TEST_PARSE_MAPPING_BASES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <assert.h>
#include <regex>
#include <string>
#include <algorithm>

using namespace std;

struct plp_elem{
    plp_elem(string cur, string ins){
        this->current_nucle = cur;
        this->insertion_nucles = ins;
    }
    string current_nucle;
    string insertion_nucles;
};

struct matrix_ele
{
    matrix_ele(int cur, int ins){
        this->current_nucle = cur;
        this->insertion_nucles = ins;
    }
    int current_nucle;
    int insertion_nucles;
};

class Matrix
{
public:
    int row_num=0, col_num=0;
    std::vector<std::vector<matrix_ele>> m;

public:
    Matrix(){};
    Matrix(int max_coverage);
    ~Matrix(){};
    void Insert(matrix_ele init_elem);
    void to_file(string filename);
};
