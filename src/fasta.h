//
// Created by Baoxing Song on 2019-04-02.
//

#ifndef WSA_FASTA_H
#define WSA_FASTA_H

#include <sstream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <iostream>
#include <string>
#include <cassert>
#include <vector>
void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences );
void charByCharCategoryRead( const std::string& filePath, std::map<std::string, std::string> & seqs, std::map<std::string, int16_t *> & weight );
#endif //WSA_FASTA_H
