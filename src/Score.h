//
// Created by Baoxing Song on 2019-04-04.
//

#ifndef WSA_SCORE_H
#define WSA_SCORE_H


#include <map>
#include <string>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <utility>

class Score {
    private:
        //std::map<int16_t, int16_t**> scoreMatrix;
        std::map<int16_t, int16_t> matchScores;
        std::map<int16_t, int16_t> mismatchScores;
        std::map<int16_t, int16_t> openPenalties;
        std::map<int16_t, int16_t> extendPenalties;
    public:
        Score(const std::string & folder);
        int16_t getScore(const int16_t & category, const char & referencer, const char & query) ;
        int16_t getOpenPenalty(const int16_t & category) ;
        int16_t getExtendPenalty(const int16_t & category) ;
};

#endif //WSA_SCORE_H
