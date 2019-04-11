/*
 * =====================================================================================
 *
 *       Filename:  fasta.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  02/04/2019 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************
implement fasta file and records operation funtion in this file



 ************************************************************************/


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
