//
// Created by Baoxing Song on 2019-04-08.
//

#include "../SuffixTree.h"
#include <stdio.h>

#include "include/gtest/gtest.h"

#include <string>
#include <iostream>
#include <sstream>

TEST(SuffixTree, c1){
    std::string seq = "AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCATGAAAGGTTATCACCCTTCTGTACTCTTATCAAGTACAACCACTCACAGAGAAGTGTGTCGTGA";
    int8_t k=7;
    int32_t seq_length = seq.size();
    std::vector<NODE> allNodes;
    CreateTree(seq.c_str(), seq_length, k, allNodes);
    std::cout << "allNodes.size() " << allNodes.size() << std::endl;
    std::cout << allNodes[0].c <<  allNodes[1].c <<  allNodes[2].c <<  allNodes[3].c <<  allNodes[4].c <<  allNodes[5].c << allNodes[6].c << allNodes[7].c << std::endl;
    std::string qSeq = "AACCCTAAACCCTAAACCCTAAACCCTAA";
    int8_t maximumMismatch=2;
    std::vector<int32_t > databasePositions;
    std::vector<int32_t > queryPositions;
    std::vector<int8_t > numMisMacth;
    findSubString( allNodes, qSeq.c_str(), qSeq.size(), k, maximumMismatch, databasePositions, queryPositions, numMisMacth );
    for( int i=0; i<databasePositions.size(); ++i ){
        std::cout << "databasePosition " << databasePositions[i] << " queryPosition " << queryPositions[i] << " numMisMacth " << (int16_t) numMisMacth[i] << std::endl;
    }
    ASSERT_EQ(0, 0);
}


TEST(SuffixTree, c2){
    std::string seq = "AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCATGAAAGGTTATCACCCTTCTGTACTCTTATCAAGTACAACCACTCACAGAGAAGTGTGTCGTGA";
    int8_t k=7;
    int32_t seq_length = seq.size();
    std::vector<NODE> allNodes;
    CreateTree(seq.c_str(), seq_length, k, allNodes);
    std::cout << "allNodes.size() " << allNodes.size() << std::endl;
    std::cout << allNodes[0].c <<  allNodes[1].c <<  allNodes[2].c <<  allNodes[3].c <<  allNodes[4].c <<  allNodes[5].c << allNodes[6].c << allNodes[7].c << std::endl;
    std::string qSeq = "AACCCTAAACCCTAAACCCTAAACCCTAA";
    int8_t maximumMismatch=2;
    std::vector<int32_t > databasePositions;
    std::vector<int32_t > queryPositions;
    std::vector<int8_t > numMisMacth;
    for( int32_t i=0; i<=(qSeq.size()-k); ++i) {
        databasePositions.clear();
        queryPositions.clear();
        numMisMacth.clear();
        findSubString(allNodes, qSeq.c_str()+i, k, maximumMismatch, databasePositions, queryPositions, numMisMacth);
        for (int j = 0; j < databasePositions.size(); ++j) {

            std::cout << "databasePosition " << databasePositions[j] << " queryPosition " << i
                      << " numMisMacth " << (int16_t)numMisMacth[j] << std::endl;
        }
    }
    ASSERT_EQ(0, 0);
}
