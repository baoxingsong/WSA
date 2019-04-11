/*
 * =====================================================================================
 *
 *       Filename:  gffToCategory.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/04/2019 19:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************
generate weigth vector by reading gff file



 ************************************************************************/


#ifndef WSA_GFFTOCATEGORY_H
#define WSA_GFFTOCATEGORY_H


#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>


void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, const std::string & outputFile);
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, std::map<std::string, int16_t *> & categories);
#endif //WSA_GFFTOCATEGORY_H
