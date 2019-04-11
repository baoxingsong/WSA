//
// Created by Baoxing Song on 2019-04-04.
//

#include "Score.h"

//
//Score::Score(const std::string & folder){
//    std::ifstream infile(folder+"/matrixList");
//    if( ! infile.good()){
//        std::cerr << "error in opening score matrix list file " << folder+"/matrixList" << std::endl;
//        exit (1);
//    }
//    std::string line="";
//    while (std::getline(infile, line)){
//        if( line.size()>0  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
//            int16_t index = std::stoi(line);
//            int16_t thisMatrix [5][5];
//            this->scoreMatrix[index] = thisMatrix;
//            std::ifstream infile2(folder+"/"+line);
//            if( ! infile2.good()){
//                std::cerr << "error in opening score matrix list file " << folder << "/" << line << std::endl;
//                exit (1);
//            }
//            std::string line2="";
//            int lineindex=0;
//            while (std::getline(infile2, line2)){
//                if( lineindex<5  ){
//                    int itermIndex=0;
//                    std::stringstream ss(line2);
//                    std::string item;
//                    while (std::getline(ss, item, '\t')) {
//                        if (item.length() > 0) {
//                            this->scoreMatrix[index][lineindex][itermIndex]=std::stoi(item);
//                        }
//                        ++itermIndex;
//                    }
//                }else if (lineindex==5){
//                    this->openPenalties[index]=std::stoi(line2);
//                }else if (lineindex==6){
//                    this->extendPenalties[index]=std::stoi(line2);
//                }
//                ++lineindex;
//            }
//            infile.close();
//        }
//    }
//    infile.close();
//}


//read score from files
Score::Score(const std::string & folder){
    std::ifstream infile(folder+"/matrixList");
    if( ! infile.good()){
        std::cerr << "error in opening score matrix list file " << folder+"/matrixList" << std::endl;
        exit (1);
    }
    std::string line="";
    while (std::getline(infile, line)){
        if( line.size()>0  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            int16_t index = std::stoi(line);
            std::ifstream infile2(folder+"/"+line);
            if( ! infile2.good()){
                std::cerr << "error in opening score matrix list file " << folder << "/" << line << std::endl;
                exit (1);
            }
            std::string line2="";
            int lineindex=0;
            while (std::getline(infile2, line2)){
                if( lineindex==0  ){
                    this->matchScores[index]=std::stoi(line2);
                }else if( lineindex==1  ){
                    this->mismatchScores[index]=std::stoi(line2);
                }else if (lineindex==2){
                    this->openPenalties[index]=std::stoi(line2);
                }else if (lineindex==3){
                    this->extendPenalties[index]=std::stoi(line2);
                }
                ++lineindex;
            }
            infile2.close();
        }
    }
    infile.close();
}


//match and mis-match query function
int16_t Score::getScore(const int16_t & category, const char & referencer, const char & query) {
    if( referencer==query ){
        return this->matchScores[category];
    }else{
        return this->mismatchScores[category];
    }
}

// open gap penalty query function
int16_t Score::getOpenPenalty(const int16_t & category) {
    return this->openPenalties[category];
}

// extend gap penalty query function
int16_t Score::getExtendPenalty(const int16_t & category) {
    return this->extendPenalties[category];
}
