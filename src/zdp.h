//
// Created by Baoxing Song on 2019-04-04.
//

#ifndef WSA_ZDP_H
#define WSA_ZDP_H

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stack>

#include "Score.h"

enum VARIANTCATEGORY {
    SNP, INSERTION, DELETION, SNPORINSERTION, SNPORDELETION, INSERTIONORDELETION, SNPORINSERTIONORDELETION
};

void zdp(const char * seqAChar, const char * seqBChar, const int32_t & lengthA, const int32_t & lengthB, const int16_t * category, std::stack<char> & A, Score & score);

#endif //WSA_ZDP_H
