#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <string>
#include "cluster.h"
#include <iomanip>
#include <random>

string read_para(char* argv[]);
void read_data(std::string filename, int nc, Cluster *clusters);
void write_data(double* score);
void file_merge(int nfiles);

#endif
