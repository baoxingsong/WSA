//
// Created by song on 9/19/18.
//


#include "include/gtest/gtest.h"
#include "../seed.h"
#include "../fasta.h"
#include "../Score.h"
#include "../gffToCategory.h"
#include <string>
#include <iostream>
#include <sstream>

TEST(seed, c1){
     std::string seqA="ATGGCGAATCTCCAGAGGGAGAGGAACCCAAGGTTGTGGCTCATAGAGAAGATATAAAGGCGCAGTTCAAGGCCGCACAACAGGTGGCGAAACCTCGTGCTCATTGGCCTCCAGAGGCCGCAAAAGCTTTCAGTGAAATATGCATTGAGAAAAGGGATTTAATCCCTTCAGGAAAAGCTAATAGTTGGGTTGGACACCCTATATCTGATGAAGTTGGATTAGAGCTTTTTCTTCTGACAAATATAGCCTACACAAGTAAACAGATGGGAAACAAGTGGTGCAACATGAGAGAATTACAACATAAAAAAGACTTGATACAAACAAGAATGTTTGATGAGCCCCACTCAGCTGCCGGGCACCTCAGACTGCCAGGACAAGAGGAAAATGAAAGTGTGAAATATCTTAAGCTAGAAGAAGAGAATGCCCGATTAAATAAACAACTGCAGGAAACGTTGACAAAAGTTAACCATCTTCAAGATCTACTTCAAAAAGCAGATGCAATACGTAGCACAATGCATAATCAAATAATGGAACTGAAGAAAAGCATACGTGTGTTTTGCCGTGTGCGACCATTGCTACAAACAGAACACAAACAGGCCAGAATCACATTTCCTGAAACATGGGGATACATTGGCCGATGTGTTTGTTTGGCACATTCAGAGTCAGAATATTCTTTCACCTGTGATAGAGTCTTTGACCATACAGCTTCACAAGAAACAATATTCAATGAAATTTCAGAATTAGTGCAGAGTGCACTTGATGGACACAAGGTGTGCATATTTGCTTATGGTCAAACCGGTTCCGGTAAGACCTACACCATGATTGGTGAGAGAGAAGGGGATGACATGGGGATAATTCCTCGCACACTTGCGAAAATATTTGAGACCATTGAAATTCGACAATCTTATGGCTCAAATCACACCGTAATGATTTCAATGTTTGAGGTTTATAATGACAAAGTTCAAGATTTATTATCCCCATTGGCACACAATGATGACTTGCAATATGTTGTTGTTAAAACAGCAACTGAGGCACTTGGATATTTGGAGAAATCATTGAAAAAGAGGTCTGTTGCCAGAACTAAACTGAATGATAGATCATCGCGAAGCCACTTTATTTTTAAGCTGTCCATATGTTCCACTGATAAGATGCTTGAAGGGGTTGTCAACATTATCGACTTAGCTGGTAGTGAACGCCCCGACAGTGATGCGTCACGTGACCTTCAAGAGGAGGCCAAGGCCATAAATTTAAGTCTACTCGATTTGGGCATTGTCCTTTCAAAAATGGAGAATGGAGAAAATAATATATCATTTAGAGGATCGACTCTGTCTCAAATTTTAGAGAGATCACTTGGAAAAGATTTTAAAACATTAATGTTGCTTAATGTTGCATCGGAAGAAAAGTTTGCTTCTGAAACATTAGCAACATTTGATTTTGCTTCGAGGATAAAATCATGCAAAAAGACAATGAGGATAATAAAAGCTAATGAACAGAAGATTGTGTTAAGAGGAAAGGAAGGATTCAAATCGCCAGGGGGTTGGCATGAAAAAGGCGGACCGTGCATAAAGTCTCCATCTATAGCAAGGAATACAAGAAAGAAAGTTATGAGCACTGCAAATAGGCAAGGAGACTATAATAGAAAGAGGAAGATGTATCCTAACGATTGA";
     std::string seqB="ATGGCGAATCTCCAGAGAGGAGAGGAACCCAAGGTTGTGGCTCATAGAGAAGATATAAAGGCGCAGTTCAAGGCCGCACAACAGGTGGCGAAACCTCGTGCTCATTGGCCTCCAGAGGCCGCAAAAGCTTTCAGTGAAATATGCATTGAGAAAAGGGATTTAATCCCTTCAGGAAAAGCTAATAGTTGGGTTGGACACCCTATATCTGATGAAGTTGGATTAGAGCTTTTTCTTCTGACAAATATAGCCTACACAAGTAAACAGATGGGAAACAAGTGGTGCAACATGAGAGAATTACAACATAAAAAAGACTTGATACAAACAAGAATGTTTGATGAGCCCCACTCAGCTGCCGGGCACCTCAGACTGCCAGGACAAGAGGAAAATGAAAGTGTGAAATATCTTAAGCTAGAAGAAGAGAATGCCCGATTAAATAAACAACTGCAGGAAACGTTGACAAAAGTTAACCATCTTCAAGATCTACTTCAAAAAGCAGATGCAATACGTAGCACAATGCATAATCAAATAATGGAACTGAAGAAAAGCATACGTGTGTTTTGCCGTGTGCGACCATTGCTACAAACAGAACACAAACAGGCCAGAATCACATTTCCTGAAACATGGGGATACATTGGCCGATGTGTTTGTTTGGCACATTCAGAGTCAGAATATTCTTTCACCTGTGATAGAGTCTTTGACCATACAGCTTCACAAGAAACAATATTCAATGAAATTTCAGAATTAGTGCAGAGTGCACTTGATGGACACAAGGTGTGCATATTTGCTTATGGTCAAACCGGTTCCGGTAAGACCTACACCATGATTGGTGAGAGAGAAGGGGATGACATGGGGATAATTCCTCGCACACTTGCGAAAATATTTGAGACCATTGAAATTCGACAATCTTATGGCTCAAATCACACCGTAATGATTTCAATGTTTGAGGTTTATAATGACAAAGTTCAAGATTTATTATCCCCATTGGCACACAATGATGACTTGCAATATGTTGTTGTTAAAACAGCAACTGAGGCACTTGGATATTTGGAGAAATCATTGAAAAAGAGGTCTGTTGCCAGAACTAAACTGAATGATAGATCATCGCGAAGCCACTTTATTTTTAAGCTGTCCATATGTTCCACTGATAAGATGCTTGAAGGGGTTGTCAACATTATCGACTTAGCTGGTAGTGAACGCCCCGACAGTGATGCGTCACGTGACCTTCAAGAGGAGGCCAAGGCCATAAATTTAAGTCTACTCGATTTGGGCATTGTCCTTTCAAAAATGGAGAATGGAGAAAATAATATATCATTTAGAGGATCGACTCTGTCTCAAATTTTAGAGAGATCACTTGGAAAAGATTTTAAAACATTAATGTTGCTTAATGTTGCATCGGAAGAAAAGTTTGCTTCTGAAACATTAGCAACATTTGATTTTGCTTCGAGGATAAAATCATGCAAAAAGACAATGAGGATAATAAAAGCTAATGAACAGAAGATTGTGTTAAGAGGAAAGGAAGGATTCAAATCGCCAGGGGGTTGGCATGAAAAAGGCGGACCGTGCATAAAGTCTCCATCTATAGCAAGGAATACAAGAAAGAAAGTTATGAGCACTGCAAATAGGCAAGGAGACTATAATAGAAAGAGGAAGATGTATCCTAACGATTGA";

//    std::string seqA="ATGGCGAATCTCCAGAGGGAGAGGAA";
//    std::string seqB="ATGGCGAATCCCAGAGGGAGAGGAA";

    /*currently this code should work for chromosome with length ~2g and does not work for longer one
     *
     * */

    std::map<std::string, std::string> b73;
    readFastaFile( "/Users/bs674/b73.fa",  b73);

    std::map<std::string, std::string> mo17;
    readFastaFile( "/Users/bs674/Mo17.fa",  mo17);

    // the initialSeedLength and initialSeedsScoreThreadsHold value was set basing the default score parameter of minimap2 match:2 mismatch:-4 opengap:-4 extendgap:-2
    // in a 8 match alignment there could be one indel in the middle, so we set initialSeedLength as 7
    // in a 7 bp alignment only one mismatch is possiable

    int lengthWantTouse=10000;
    int16_t category[10000] ={ 0 };
    int16_t initialSeedLength=7;
    int16_t initialSeedsScoreThreadsHold=6;

    Score score("/Users/bs674/scoreMatrix");
    std::string outputFile="/Users/bs674/WSA_align_seed_c1";
    int64_t numberOfSeedsForChain = 10000;
    seq2seed ( b73["1"].substr(0, lengthWantTouse).c_str(),  mo17["1"].substr(0, lengthWantTouse).c_str(), lengthWantTouse, lengthWantTouse, category, initialSeedLength,
            initialSeedsScoreThreadsHold, outputFile, 12, score, numberOfSeedsForChain);

//    const char * seqAChar = seqA.c_str();
//    const char * seqBChar = seqB.c_str();
//    zdp(seqAChar, seqBChar, 10, 15, category, score);

    ASSERT_EQ(0, 0);
}
TEST(sore, c1){
    Score score("/Users/bs674/scoreMatrix");
    std::vector<int16_t> cs;
    cs.push_back(0);
    cs.push_back(1);
    cs.push_back(10);
    cs.push_back(70);
    cs.push_back(85);
    cs.push_back(100);
    for( int16_t  c : cs ){
        std::cout << score.getScore(c, 'A', 'A') << " " << score.getScore(c, 'A', 'C') << " " << score.getOpenPenalty(c) << " " << score.getExtendPenalty(c) << std::endl;
    }
}
TEST(gffToCategory, c1){
    std::map<std::string, std::string> b73;
    readFastaFile( "/Users/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",  b73);
    std::map<std::string, int32_t> chrSize;
    for( std::map<std::string, std::string>::iterator it=b73.begin(); it != b73.end(); ++it){
        chrSize[it->first] = it->second.size();
    }
    std::string gffFile = "/Users/bs674/Zea_mays.B73_RefGen_v4.42.gff3";
    std::string outputFile = "/Users/bs674/category";
    readGffFileWithEveryThing (gffFile,  chrSize, outputFile);
}

TEST(seed, c2){
    std::map<std::string, std::string> b73;
    //readFastaFile( "/Users/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",  b73);
    readFastaFile( "/Users/bs674/Zea_mays.AGPv3.31.dna.genome.fa",  b73);

    std::map<std::string, std::string> mo17;
    readFastaFile( "/Users/bs674/Mo17.fa",  mo17);
    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = b73.begin(); it!=b73.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    //readGffFileWithEveryThing ( "/Users/bs674/Zea_mays.B73_RefGen_v4.42.gff3",  chrSize, weight);
    readGffFileWithEveryThing ( "/Users/bs674/Zea_mays.AGPv3.31.gff3",  chrSize, weight);
    int lengthWantTouse=100000;

    int8_t initialSeedLength=7;
    int8_t initialSeedsScoreThreadsHold=6;

    Score score("/Users/bs674/scoreMatrix");
    std::cout << "start alignment" << std::endl;
    std::string outputFile="/Users/bs674/WSA_align_seed_c3";
    int64_t numberOfSeedsForChain = 0;
    seq2seed ( b73["1"].c_str()+4000,  mo17["1"].c_str()+624000, lengthWantTouse, lengthWantTouse, weight["1"]+4000,
               initialSeedLength, initialSeedsScoreThreadsHold, outputFile, 50, score, numberOfSeedsForChain);

    for( std::map<std::string, int16_t *>::iterator it=weight.begin(); it != weight.end(); ++it ){
        delete[] it->second;
    }
    ASSERT_EQ(0, 0);
}

TEST(seed, c3){ // this one does not work so well, since my code could not deal with very long k-mersss
    std::map<std::string, std::string> b73;
    //readFastaFile( "/Users/bs674/Zea_mays.B73_RefGen_v4.dna.toplevel.fa",  b73);
    readFastaFile( "/Users/bs674/Zea_mays.AGPv3.31.dna.genome.fa",  b73);

    std::map<std::string, std::string> mo17;
    readFastaFile( "/Users/bs674/Mo17.fa",  mo17);
    std::map<std::string, int32_t>  chrSize;
    for( std::map<std::string, std::string>::iterator it = b73.begin(); it!=b73.end(); ++it ){
        chrSize[it->first] = it->second.size();
    }
    std::map<std::string, int16_t *> weight;
    //readGffFileWithEveryThing ( "/Users/bs674/Zea_mays.B73_RefGen_v4.42.gff3",  chrSize, weight);
    readGffFileWithEveryThing ( "/Users/bs674/Zea_mays.AGPv3.31.gff3",  chrSize, weight);
    //int lengthWantTouse=300000;

    int8_t initialSeedLength=8;
    int8_t initialSeedsScoreThreadsHold=8;

    Score score("/Users/bs674/scoreMatrix");
    std::cout << "start alignment" << std::endl;
    std::string outputFile="/Users/bs674/WSA_align_seed_c3";
    int64_t numberOfSeedsForChain = 10000;

    seq2seed ( b73["1"].c_str(),  mo17["1"].c_str(),
               b73["1"].size(), mo17["1"].size(), weight["1"],
            initialSeedLength, initialSeedsScoreThreadsHold, outputFile, 2000, score, numberOfSeedsForChain);

    for( std::map<std::string, int16_t *>::iterator it=weight.begin(); it != weight.end(); ++it ){
        delete[] it->second;
    }
    ASSERT_EQ(0, 0);
}
