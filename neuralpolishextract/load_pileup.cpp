//
// Created by lab on 2020/7/17.
//

#include "load_pileup.h"

Matrix::Matrix(int max_coverage) {
    row_num = max_coverage;
    m.resize(row_num);
}

void Matrix::Insert(matrix_ele init_elem) {
    col_num++;
    for (int i = 0; i < row_num; i++) {
        m[i].push_back(init_elem);
    }
}

void Matrix::to_file(string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i < row_num; i++) {
        for (vector<matrix_ele>::iterator iter = m[i].begin(); iter != m[i].end(); iter++) {
            outfile << iter->current_nucle << "," << iter->insertion_nucles;
            if (iter != m[i].end() - 1)
                outfile << ",";
            else
                outfile << "\n";
        }
    }
    outfile.close();
}

void SplitString(string &s, vector<string> &v, string &c) {
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while (string::npos != pos2) {
        v.push_back(s.substr(pos1, pos2 - pos1));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1));
}

void parse_mapping_bases(string &mapping_bases, regex &r, vector<plp_elem> &parsed_vector, string &ref_base, int max_coverage) {
    sregex_iterator it(mapping_bases.begin(), mapping_bases.end(), r);
    sregex_iterator end;
    int parsed_count = 0;
    while (it != end) {
        if(parsed_count>max_coverage) break;
        string current_nucle, insertion_nucles;
        if (it->str()[0] == '^' || it->str()[0] == '$') {
            it++;
            continue;
        } else if (it->str() == "." || it->str() == ",") {
            // 匹配
            current_nucle = ref_base;
            insertion_nucles = "B";
            plp_elem t(current_nucle, insertion_nucles);
            parsed_vector.push_back(t);
        } else if (it->str() == "*") {
            // 缺失
            current_nucle = "D";
            insertion_nucles = "B";
            plp_elem t(current_nucle, insertion_nucles);
            parsed_vector.push_back(t);
        } else if (it->str()[0] == '+') {
            // 插入
            int indel_len = stoi(it->str(2));
            string indel_nucles = it->str(3);
            if (indel_nucles.size() == indel_len) {
                replace(indel_nucles.begin(), indel_nucles.end(), 'N', '\0');
                indel_nucles = indel_nucles.size() > 5 ? indel_nucles.substr(0, 5) : indel_nucles;
                for_each(indel_nucles.begin(), indel_nucles.end(), [](char &c) { c = toupper(c); });
                parsed_vector.back().insertion_nucles = indel_nucles;
            } else {
                string true_insertion, comming_mismatches;
                true_insertion = indel_nucles.substr(0, indel_len);
                comming_mismatches = indel_nucles.substr(indel_len, indel_nucles.size() - indel_len);

                replace(true_insertion.begin(), true_insertion.end(), 'N', '\0');
                true_insertion = true_insertion.size() > 5 ? true_insertion.substr(0, 5) : true_insertion;
                for_each(true_insertion.begin(), true_insertion.end(), [](char &c) { c = toupper(c); });
                parsed_vector.back().insertion_nucles = true_insertion;

                for (int mi = 0; mi < comming_mismatches.size(); mi++) {
                    if (comming_mismatches[mi] == 'N' || comming_mismatches[mi] == 'n') {
                        current_nucle = ref_base;
                        insertion_nucles = "B";
                        plp_elem t(current_nucle, insertion_nucles);
                        parsed_vector.push_back(t);
                    } else {
                        current_nucle = toupper(comming_mismatches[mi]);
                        insertion_nucles = "B";
                        plp_elem t(current_nucle, insertion_nucles);
                        parsed_vector.push_back(t);
                    }
                }
            }
        } else if (it->str()[0] == '-') {
            // 缺失
            int indel_len = stoi(it->str(5));
            string indel_nucles = it->str(6);
            if (indel_nucles.size() != indel_len) {
                string comming_mismatches;
                comming_mismatches = indel_nucles.substr(indel_len, indel_nucles.size() - indel_len);
                for (int mi = 0; mi < comming_mismatches.size(); mi++) {
                    if (comming_mismatches[mi] == 'N' || comming_mismatches[mi] == 'n') {
                        current_nucle = ref_base;
                        insertion_nucles = "B";
                        plp_elem t(current_nucle, insertion_nucles);
                        parsed_vector.push_back(t);
                    } else {
                        current_nucle = toupper(comming_mismatches[mi]);
                        insertion_nucles = "B";
                        plp_elem t(current_nucle, insertion_nucles);
                        parsed_vector.push_back(t);
                    }
                }

            }
        } else {
            assert(it->str().size() == 1);
            // 不匹配
            if (it->str() == "N" || it->str() == "n") {
                current_nucle = ref_base;
                insertion_nucles = "B";
                plp_elem t(current_nucle, insertion_nucles);
                parsed_vector.push_back(t);
            } else {
                current_nucle = toupper(it->str()[0]);
                insertion_nucles = "B";
                plp_elem t(current_nucle, insertion_nucles);
                parsed_vector.push_back(t);
            }
        }
        ++it;
        parsed_count++;
    }
}

void fill_in_matrix(vector<plp_elem> &align_base_vector, vector<string> &mapping_readnames, Matrix &mat,
                    map<string, int> &encoder, map<string, int> &read_name_idx, int &row_count) {
    string read_name;
    int row_idx;
    for (int i = 0; i < mapping_readnames.size(); i++) {
        read_name = mapping_readnames[i];
        if (read_name_idx.find(read_name) != read_name_idx.end()) {
            //存在
            row_idx = read_name_idx.at(read_name);
        } else {
            //不存在
            row_idx = row_count;
            read_name_idx.insert(pair<string, int>(read_name, row_count++));
        }
        if (row_idx >= mat.row_num) continue;
        mat.m[row_idx].back().current_nucle = encoder.at(align_base_vector[i].current_nucle);
        if (align_base_vector[i].insertion_nucles != "B")
            mat.m[row_idx].back().insertion_nucles = encoder.at(align_base_vector[i].insertion_nucles);

    }
}

pair<vector<string>,vector<Matrix>> read_pileup(string filename, uint32_t window_length, uint32_t mat_coverage, map<string, int> &encoder) {
    ifstream fin;
    string line;
    string contig_name, ref_base, mapping_bases, mapping_quals, mapping_readnames;
    int ref_pos, coverage, pre_pos = -1;
    string readname_sep = ",";
    string pattern = "(\\+)([0123456789]+)([ACGTNacgtn]+)|(\\-)([0123456789]+)([ACGTNacgtn]+)|\\^.|\\$|\\.|,|[ACGTN]|[acgtn]|\\*";
    regex r_pattern(pattern);
    int max_coverage = mat_coverage;
    std::string blank_base = "B";
    matrix_ele init_elem(encoder.at(blank_base), encoder.at(blank_base));
    map<string, int> read_name_idx;
    int row_count=0;
    int region_start, region_end;
    Matrix mat(max_coverage);
    vector<Matrix> mat_vector;
    vector<string> region_str_vector;
    pair<vector<string>,vector<Matrix>> result;

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
        contig_name = fields[0];
        ref_pos = stoi(fields[1]);
        ref_base = fields[2];
        coverage = stoi(fields[3]);
        mapping_bases = fields[4];
        mapping_quals = fields[5];
        mapping_readnames = fields[6];
        vector<plp_elem> align_base_vector;
        parse_mapping_bases(mapping_bases, r_pattern, align_base_vector, ref_base, max_coverage);
        vector<string> mapping_readnames_vector,used_readnames_vector;
        SplitString(mapping_readnames, mapping_readnames_vector, readname_sep);
        //检查mapping reads和reads name数量是否一样
        for(int vi=0;vi<align_base_vector.size();vi++){
            used_readnames_vector.push_back(mapping_readnames_vector[vi]);
        }
        assert(used_readnames_vector.size() == align_base_vector.size());
        if (pre_pos == -1) {//处理文件第一行
            pre_pos = ref_pos - 1;
            region_end = region_start = ref_pos;
        }
        if (mat.col_num == window_length || ref_pos != pre_pos + 1) {
            if (ref_pos != pre_pos + 1) cout << "position: " << ref_pos << " is disconnected." << endl;
            // 已处理64列，或者遇到断开位点
            string region_str = contig_name + ":" + to_string(region_start) + "-" + to_string(region_end);
            mat.row_num = row_count < max_coverage ? row_count : max_coverage;//只使用实际层数的数据
//            mat.to_file(region_str);
            mat_vector.push_back(mat);
            region_str_vector.push_back(region_str);
            mat = Matrix(max_coverage);
            row_count = 0;
            read_name_idx.clear();
            region_end = region_start = ref_pos;
        }
        mat.Insert(init_elem);
        region_end = ref_pos;
        fill_in_matrix(align_base_vector, used_readnames_vector, mat, encoder, read_name_idx, row_count);
        pre_pos = ref_pos;
    }
    if(mat.col_num>0){
        string region_str = contig_name + ":" + to_string(region_start) + "-" + to_string(region_end);
        mat.row_num = row_count < max_coverage ? row_count : max_coverage;//只使用实际层数的数据
//            mat.to_file(region_str);
        mat_vector.push_back(mat);
        region_str_vector.push_back(region_str);
    }
    result.first = region_str_vector;
    result.second = mat_vector;
    return result;
}