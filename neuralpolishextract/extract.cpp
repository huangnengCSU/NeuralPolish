//
// Created by lab on 2020/7/17.
//

#include "extract.h"

using namespace std;

const char *encoder_file = "encoder.txt";
const char *decoder_file = "decoder.txt";

bool get_all_files(const string &dir_in, vector<string> &files) {
    if (dir_in.empty()) {
        return false;
    }
    struct stat s;
    stat(dir_in.c_str(), &s);
    if (!S_ISDIR(s.st_mode)) {
        return false;
    }
    DIR *open_dir = opendir(dir_in.c_str());
    if (NULL == open_dir) {
        std::exit(EXIT_FAILURE);
    }
    dirent *p = nullptr;
    while ((p = readdir(open_dir)) != nullptr) {
        struct stat st;
        if (p->d_name[0] != '.') {
            //因为是使用devC++ 获取windows下的文件，所以使用了 "\" ,linux下要换成"/"
            std::string name = dir_in + std::string("/") + std::string(p->d_name);
            stat(name.c_str(), &st);
            if (S_ISDIR(st.st_mode)) {
                get_all_files(name, files);
            } else if (S_ISREG(st.st_mode)) {
                files.push_back(name);
            }
        }
    }
    closedir(open_dir);
    return true;
}

std::map<std::string, int> load_encoder(const char *encoder_file) {
    std::map<std::string, int> encoder;
    string kmer;
    int value;
    ifstream infile;
    infile.open(encoder_file);
    while (!infile.eof()) {
        infile >> kmer >> value;
        encoder.insert(std::pair<std::string, int>(kmer, value));
    }
    return encoder;
}

std::map<int, std::string> load_decoder(const char *decoder_file) {
    std::map<int, std::string> decoder;
    string kmer;
    int value;
    ifstream infile;
    infile.open(decoder_file);
    while (!infile.eof()) {
        infile >> value >> kmer;
        decoder.insert(std::pair<int, std::string>(value, kmer));
    }
    return decoder;
}

struct CmdArgs {
    char *pile_dir = nullptr;
    char *out_dir = nullptr;
    uint32_t window_size = 64;
    uint32_t num_threads = 1;
    uint32_t max_coverage = 40;
};

CmdArgs getOpt(int argc, char *argv[]) {
    int opt;
    const char *optstring = "p:w:t:o:c:";
    CmdArgs args;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        if ((char) opt == 'p')
            args.pile_dir = optarg;
        else if ((char) opt == 'o')
            args.out_dir = optarg;
        else if ((char) opt == 'w') {
            string w_size = optarg;
            args.window_size = atoi(w_size.c_str());
        } else if ((char) opt == 't') {
            string n_threads = optarg;
            args.num_threads = atoi(n_threads.c_str());
        } else if ((char) opt == 'c') {
            string coverage = optarg;
            args.max_coverage = atoi(coverage.c_str());
        }
    }
    return args;
}

int main(int argc, char *argv[]) {
    CmdArgs args;
    args = getOpt(argc, argv);
    int isCreate = mkdir(args.out_dir, S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
    std::map<std::string, int> encoder = load_encoder(encoder_file);
    vector<string> filelist;
    string dir_in = args.pile_dir;
    string dir_out = args.out_dir;
    get_all_files(dir_in, filelist);
    generate_features(filelist, args.window_size, encoder, args.num_threads, args.max_coverage, dir_out);
    return 0;
}