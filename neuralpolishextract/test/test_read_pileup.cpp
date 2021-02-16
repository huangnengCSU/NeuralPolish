//
// Created by lab on 2020/7/17.
//
#include "test_read_pileup.h"
void test_read_pileup(string filename) {
    ifstream fin;
    string line;
    string contig_name, mapping_bases, mapping_readnames;
    int ref_pos, coverage;
    char ref_base;
    fin.open(filename);
    if (!fin.is_open()) {
        cout << "未成功打开文件" << endl;
    }
    while (getline(fin, line)) {
        stringstream ss(line);
        string tmp;
        vector<string> fields;
        while (getline(ss, tmp, '\t')) {
            fields.push_back(tmp);
        }
        for (int i = 0; i < 4; i++) {
            cout << fields[i] << "\t";
        }
        cout << endl;
    }
}

int main() {
    string filename = "/home/user/code/neuralpolishextract/ecoli_pileups/contig_1:3440001-3450000.mpileup";
    test_read_pileup(filename);
    return 0;
}