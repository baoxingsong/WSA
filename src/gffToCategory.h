//
// Created by Baoxing Song on 2019-04-06.
//

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
