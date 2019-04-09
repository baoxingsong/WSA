//
// Created by Baoxing Song on 2019-04-02.
//

#ifndef WSA_SEED_H
#define WSA_SEED_H
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include "Score.h"
#include "zdp.h"

void seq2seed ( const std::string & seqA, const std::string & seqB, const int16_t * category,
                const int8_t & initialSeedLength, const int8_t & initialSeedsScoreThreadsHold,
                const std::string & outputFile, int8_t miniseedLength, Score & score,
                const int64_t & numberOfSeedsForChain);

#endif //WSA_SEED_H
