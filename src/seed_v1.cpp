//
// Created by Baoxing Song on 2019-04-02.
//

#include <cstdint>
#include "seed.h"

/*

struct MergedSeed {
    int32_t length; // this is only used for the InitialSeed extend result
    int32_t positionAStart;
    int32_t positionBStart;
    int64_t score;

    int32_t positionAEnd;
    int32_t positionBEnd;
    std::string align; // if align is "", it means every should be M, the base align is indicated using the cigar letters M D I
    bool everUsed;
};


//O2
void creatSeed(const char * a, const char * b, const int8_t & seedLength, int8_t & score){ // never use initial seed length longer than 126
    if( *(a) == *(b) && *(a+seedLength-1) == *(b+seedLength-1) ){
        score=2;
    }else{
        score=0;
        return;
    }
    for ( int8_t i=1; i < seedLength-1; ++i ) {
        score += (*(a+i) == *(b+i));
    }
}

//computational complexity should be o1
void forwardExtendSeed(const char * a, const char * b, const int32_t & sequenceLengthA, const int32_t & sequenceLengthB, const int8_t & seedLength, const int8_t & seedScore, int32_t & extendedSeedLength, int32_t & extendSeedScore){
    extendedSeedLength = seedLength;
    extendSeedScore = seedScore;
    for ( int32_t i=seedLength; i < sequenceLengthA && i < sequenceLengthB; ++i ) {
        if ( *(a+i) == *(b+i) ){
            ++extendSeedScore;
            ++extendedSeedLength;
        }else{
            return;
        }
    }
}

void megerSeedInsertion ( const int & seedScore1, const int & seedScore2, const int32_t & position1AEnd,
                 const int32_t & position2AStart, const int32_t & position1BEnd, const int32_t & position2BStart,
                 const std::string & align1, const std::string & align2, const int16_t * category,
                 Score & score, int8_t & couldMerge, int32_t & newScore,
                 std::string & align ) {
   couldMerge=0;
   if( (position1BEnd+1) > position2BStart ){
       return;
   } else if (  (position1BEnd+1) < position2BStart ){
        if( (position2BStart-position1BEnd) > seedScore1  ){
            couldMerge=2;
            return;
        }
        if(  (position2BStart-position1BEnd) > seedScore2 ){
            return;
        }
        int16_t openGapPenalty, extendGapPenalty; // here using the smaller indel penalty
        if( score.getOpenPenalty(*(category+position1AEnd)) > score.getOpenPenalty(*(category+position2AStart)) ){
            openGapPenalty=score.getOpenPenalty(*(category+position2AStart));
            extendGapPenalty=score.getExtendPenalty(*(category+position2AStart));
        }else{
            openGapPenalty=score.getOpenPenalty(*(category+position1AEnd));
            extendGapPenalty=score.getExtendPenalty(*(category+position1AEnd));
        }
        if( openGapPenalty > seedScore1 || openGapPenalty > seedScore2 ){
            return;
        }
        int32_t thisScore = seedScore1 + seedScore2 - openGapPenalty - (position2BStart-(position1BEnd+1))*extendGapPenalty;  // assume we have a seed 010-010 the first seed score is 2 and the second seed score is 2 and the merged seed score is 2 (0 means match, 1 means mis-match and - means indel)
        if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
            newScore = thisScore;
            couldMerge = 1;
        }else{
            return;
        }
        align = align1 + std::string((position2BStart-(position1BEnd+1)), 'I') + align2;
    }
}

void megerSeedDeletion ( const int & seedScore1, const int & seedScore2, const int32_t & position1AEnd,
                 const int32_t & position2AStart, const int32_t & position1BEnd, const int32_t & position2BStart,
                 const std::string & align1, const std::string & align2, const int16_t * category,
                 Score & score, int8_t & couldMerge, int32_t & newScore,
                 std::string & align ) {
    couldMerge = 0;
    if( position1AEnd+1 > position2AStart ){
        return;
    }else if( ( position1AEnd+1 < position2AStart) ){

        if( (position2AStart- position1AEnd) > seedScore1 ){ // here position2AStart- position1AEnd is gapsize + 1, since we have both OpenPenalty and ExtendPenalty, so it should work
            couldMerge=2;
            return;
        }
        if( (position2AStart- position1AEnd) > seedScore2 ){ // here position2AStart- position1AEnd is gapsize + 1, since we have both OpenPenalty and ExtendPenalty, so it should work
            return;
        }
        int32_t gapPenalty; // here using the smaller indel penalty
        int32_t miniScore = seedScore1 < seedScore2 ? seedScore1 : seedScore2;
        if( score.getOpenPenalty(*(category+position1AEnd+1)) > score.getOpenPenalty(*(category+position2AStart-1)) ){
            gapPenalty=score.getOpenPenalty(*(category+position2AStart-1));
            for( int i=0; position2AStart-1-i > position1AEnd; ++i){
                gapPenalty += score.getExtendPenalty(*(category+position2AStart-1-i));
                if( gapPenalty > miniScore ){ // if find something not worth going one, and stop it could make it faster
                    return;
                }
            }
        }else{
            gapPenalty=score.getOpenPenalty(*(category+position1AEnd+1));
            for( int i=0; position1AEnd+1+i< position2AStart; ++i){
                gapPenalty += score.getExtendPenalty(*(category+position1AEnd+1+i));
                if( gapPenalty > miniScore ){
                    return;
                }
            }
        }
        int32_t thisScore = seedScore1 + seedScore2 - gapPenalty;
        //int32_t thisScore = seedScore1 + seedScore2 - position2AStart + position1AEnd;
        if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
            newScore = thisScore;
            couldMerge = 1;
        }else{
            return;
        }
        align = align1+ std::string((position2AStart-(position1AEnd+1)), 'D') +align2;
    }
}

// the computational complexity of this step should be o1
// only apply this function on seed without indel
void removeExtendSeedsDuplications(std::vector<MergedSeed> & allExtendSeeds){
    std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](MergedSeed a, MergedSeed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    int i, j;
    for ( i=0, j=0; i< allExtendSeeds.size(); ++i ){
        if( ( (allExtendSeeds[j].positionAStart+allExtendSeeds[j].length) >= (allExtendSeeds[i].positionAStart+allExtendSeeds[i].length) )
         && allExtendSeeds[j].positionBStart <= allExtendSeeds[i].positionBStart && ( (allExtendSeeds[j].positionBStart+allExtendSeeds[j].length) >= (allExtendSeeds[i].positionBStart+allExtendSeeds[i].length )) ){

        }else{
            ++j;
            if( i != j ){
                allExtendSeeds[j]=allExtendSeeds[i];
            }
        }
    }
    ++j;
    allExtendSeeds.resize(j);
}


//o2
void longestPath_V1 (std::vector<MergedSeed> & allMergedSeeds, std::vector<MergedSeed> & chain){
    std::sort(allMergedSeeds.begin(), allMergedSeeds.end(), [](MergedSeed a, MergedSeed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    std::cout << "begin to chain" << std::endl;
    int64_t maxSore = 0;
    int32_t bestEnd = 0;
//    std::cout << allMergedSeeds[0].score << " " << allMergedSeeds[0].positionAStart << " " << allMergedSeeds[0].positionAEnd << " " << allMergedSeeds[0].positionBStart << " " << allMergedSeeds[0].positionBEnd << std::endl;
    std::vector<int64_t> scoreArray (allMergedSeeds.size()); // arrays of scores
    std::vector<int32_t> prev (allMergedSeeds.size());  // index of previous node in longest path
    scoreArray[0] = allMergedSeeds[0].score;
    prev[0] = -1;

    if (scoreArray[0] > maxSore){
        bestEnd = 0;
        maxSore = scoreArray[0];
    }
    for (int32_t idx = 1; idx < allMergedSeeds.size(); ++idx) {
        scoreArray[idx] = allMergedSeeds[idx].score;
        prev[idx] = -1;
        for (int32_t jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
            // Because we swapped asm/query start position so that inversions were all increasing,
            // we should always be on the diagonal.  If not, then we filter it.
            // This gets rid of the noise, while preserving the inversions on
            // the diagonal
            // Are only looking at positions previous to our current "idx" position
            if (allMergedSeeds[jdx].positionAEnd < allMergedSeeds[idx].positionAStart &&
                allMergedSeeds[jdx].positionBEnd < allMergedSeeds[idx].positionBStart) {
                if ((scoreArray[jdx] + allMergedSeeds[idx].score) > scoreArray[idx]) {
                    scoreArray[idx] = scoreArray[jdx] + allMergedSeeds[idx].score;
                    prev[idx] = jdx;
                }
            }
        }
        if (scoreArray[idx] > maxSore){
            bestEnd = idx;
            maxSore = scoreArray[idx];
        }
    }
    int idx=bestEnd; // bestEnd is where to stop the longest path
    chain.push_back(allMergedSeeds[idx]);
    int jdx = prev[idx]; // prev[] is index on the longest path
    while( jdx>=0 ){
        chain.push_back(allMergedSeeds[jdx]);
        jdx=prev[jdx];
    }// Reversing the order
    std::reverse(std::begin(chain), std::end(chain));
}


void longestPath (std::vector<MergedSeed> & allMergedSeeds, std::vector<MergedSeed> & chain){
    // there is problem with the idea behind this function
    // since the Bstart and Bend could not be well orginazed ofr the used set
    std::sort(allMergedSeeds.begin(), allMergedSeeds.end(), [](MergedSeed a, MergedSeed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    int64_t maxSore = 0;
    int32_t bestEnd = 0;
    std::vector<int64_t> scoreArray (allMergedSeeds.size()); // arrays of scores
    std::vector<int32_t> prev (allMergedSeeds.size());  // index of previous node in longest path
    scoreArray[0] = allMergedSeeds[0].score;
    prev[0] = -1;

    if (scoreArray[0] > maxSore){
        bestEnd = 0;
        maxSore = scoreArray[0];
    }
    std::set<int32_t> used;
    used.insert(0);
    int32_t toRemoveFromUsed;
    for (int32_t idx = 1; idx < allMergedSeeds.size(); ++idx) {
        scoreArray[idx] = allMergedSeeds[idx].score;
        prev[idx] = -1;
        toRemoveFromUsed = -1;
        for (int32_t jdx : used) {
            if (allMergedSeeds[jdx].positionAEnd < allMergedSeeds[idx].positionAStart &&
                allMergedSeeds[jdx].positionBEnd < allMergedSeeds[idx].positionBStart) {
                if ((scoreArray[jdx] + allMergedSeeds[idx].score) > scoreArray[idx]) {
                    scoreArray[idx] = scoreArray[jdx] + allMergedSeeds[idx].score;
                    prev[idx] = jdx;
                    toRemoveFromUsed = prev[jdx];
                }
            }
        }
        if( toRemoveFromUsed > -1 ){
            used.erase(toRemoveFromUsed);
        }
        used.insert(idx);
        if (scoreArray[idx] > maxSore){
            bestEnd = idx;
            maxSore = scoreArray[idx];
        }
    }

    int idx=bestEnd; // bestEnd is where to stop the longest path
    chain.push_back(allMergedSeeds[idx]);
    int jdx = prev[idx]; // prev[] is index on the longest path
    while( jdx>=0 ){
        chain.push_back(allMergedSeeds[jdx]);
        jdx=prev[jdx];
    }// Reversing the order
    std::reverse(std::begin(chain), std::end(chain));
}


void outputAlignment( const char * a, const char * b, const int & startA, const int & startB, const std::string & align ){
    std::string aa;
    std::string ab;
    int ai=startA;
    int bi=startB;
    for ( int i=0; i<align.size(); ++i ){
        if( align[i] == 'M' ){
            aa += *(a+ai);
            ab += *(b+bi);
            ++ai;
            ++bi;
        }else if( align[i] == 'D' ){
            aa += *(a+ai);
            ab += '-';
            ++ai;
        }else if( align[i] == 'I' ){
            aa += '-';
            ab += *(b+bi);
            ++bi;
        }
    }
    std::cout << aa << std::endl << ab << std::endl;
}
void outputAlignment( const char * a, const char * b, const int & startA, const int & startB, const int32_t & length ){
    std::string align="";
    for( int i =0; i<length; ++i ){
        align += "M";
    }
    outputAlignment(  a, b, startA, startB, align );
}

void seq2seed ( const std::string & seqA, const std::string & seqB, const int16_t * category, const int8_t & initialSeedLength, const int8_t & initialSeedsScoreThreadsHold, int8_t miniseedLength, Score & score){
    miniseedLength--;
    std::cout << seqA << std::endl << seqB << std::endl;

    int32_t lengthSeqA = seqA.size();
    int32_t lengthSeqB = seqB.size();
    const char * seqAChar = seqA.c_str();
    const char * seqBChar = seqB.c_str();
    int8_t initialSeedsScore;
    int32_t i, j;
    int32_t extendedSeedLength;
    int32_t extendSeedScore;
    std::vector<MergedSeed> allExtendSeeds;
    for( i=0; i<=(lengthSeqA-initialSeedLength); ++i ){
        for( j=0; j<=(lengthSeqB-initialSeedLength); ++j ){
            creatSeed(seqAChar+i, seqBChar+j, initialSeedLength, initialSeedsScore);
            if( initialSeedsScore >= initialSeedsScoreThreadsHold ){
                forwardExtendSeed( seqAChar+i, seqBChar+j, lengthSeqA-i,
                                   lengthSeqB-j, initialSeedLength, initialSeedsScore, extendedSeedLength, extendSeedScore);
                MergedSeed extendSeed={extendedSeedLength, i, j, extendSeedScore};
                allExtendSeeds.push_back(extendSeed);
                j += initialSeedLength-1;
            }
        }
    }
    std::cout << "Initial Seeds " << allExtendSeeds.size() << std::endl;

    removeExtendSeedsDuplications(allExtendSeeds); // remove duplication records
    std::cout << "forward extend of seeds remove duplication " << allExtendSeeds.size() << std::endl;

    // update seed score according to their category
    int64_t newScore;
    for( i=0; i<allExtendSeeds.size(); ++i ){
        newScore=0;
        for( j=0; j< allExtendSeeds[i].length; ++j){
            newScore += score.getScore( *(category+j+allExtendSeeds[i].positionAStart), *(seqAChar+j+allExtendSeeds[i].positionAStart), *(seqBChar+j+allExtendSeeds[i].positionBStart) ); // the category share the same coordinate with reference sequence
        }
        allExtendSeeds[i].score=newScore;
    }
    time_t my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update seed score according to their category " << allExtendSeeds.size() << std::endl;


    // the following code could be dangerous, should be really careful
    // since there are INDELs, so should not remove duplication by comparing the start and end position without consider the alignment
    //std::vector<MergedSeed> allMergedSeeds(allExtendSeeds.size());

    for( i=0; i<allExtendSeeds.size(); ++i ){
        allExtendSeeds[i].positionAEnd = allExtendSeeds[i].positionAStart+allExtendSeeds[i].length-1;
        allExtendSeeds[i].positionBEnd = allExtendSeeds[i].positionBStart+allExtendSeeds[i].length-1;
        allExtendSeeds[i].align=std::string(allExtendSeeds[i].length, 'M');
        allExtendSeeds[i].everUsed=false;
    }
    std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](MergedSeed a, MergedSeed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "extend transformed to megered seeds " << allExtendSeeds.size() << std::endl;


    // prepare special data structure for seeding up purpose
    std::vector<std::vector<MergedSeed *>> referencePositionAsIndex(lengthSeqA+1); // plus for end based entry
    size_t allMergedSeedsSize = allExtendSeeds.size();
    for( i=0; i<allMergedSeedsSize; ++i ){
        referencePositionAsIndex[allExtendSeeds[i].positionAStart].push_back(&allExtendSeeds[i]);
    }

    int32_t megerScore;
    int8_t couldMerge=2;
    std::string align;
    for( i=0; i<allMergedSeedsSize; ++i ){
        if( !allExtendSeeds[i].everUsed ){
            referencePositionAsIndex:
            if( referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1].size()>0 ){
                for (j=0; j < referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1].size(); ++j) {
                    if( !referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->everUsed ){
                        megerSeedInsertion(allExtendSeeds[i].score, referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->score, allExtendSeeds[i].positionAEnd,
                                  referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->positionAStart, allExtendSeeds[i].positionBEnd,
                                  referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->positionBStart,
                                           allExtendSeeds[i].align, referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->align, category,
                                  score, couldMerge, megerScore, align);
                        if( 2 == couldMerge ){
                            break;
                        } else if (couldMerge == 1) {
                            allExtendSeeds[i].positionBEnd = referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->positionBEnd;
                            referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->everUsed=true;
                            allExtendSeeds[i].positionAEnd = referencePositionAsIndex[allExtendSeeds[i].positionAEnd+1][j]->positionAEnd;
                            allExtendSeeds[i].align = align;
                            allExtendSeeds[i].score = megerScore;
                            goto referencePositionAsIndex;
                        }
                    }
                }
            }
        }
    }
    referencePositionAsIndex.clear();
    for( i=j=0; i<allExtendSeeds.size(); ++i ){
        if( !allExtendSeeds[i].everUsed ){ // should be kept
            if( i != j ){
                allExtendSeeds[j]=allExtendSeeds[i];
            }
            ++j;
        }
    }
    std::cout << "referencePositionAsIndex shrink allMergedSeeds " << allExtendSeeds.size() << std::endl;
    allMergedSeedsSize=allExtendSeeds.size();
    std::vector<std::vector<MergedSeed *>> queryPositionAsIndex(lengthSeqB+1);
    for( i=0; i<allMergedSeedsSize; ++i ){
        queryPositionAsIndex[allExtendSeeds[i].positionBStart].push_back(&allExtendSeeds[i]);
    }
    for( i=0; i<allMergedSeedsSize; ++i ){
        if( !allExtendSeeds[i].everUsed ){
            queryPositionAsIndex:
            if( queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1].size() > 0 ){
                for (j=0; j < queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1].size(); ++j) {
                    if( !queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1][j]->everUsed ) {
                        megerSeedDeletion(allExtendSeeds[i].score, queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1][j]->score, allExtendSeeds[i].positionAEnd,
                                          queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1][j]->positionAStart, allExtendSeeds[i].positionBEnd,
                                          queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1][j]->positionBStart,
                                          allExtendSeeds[i].align, queryPositionAsIndex[allExtendSeeds[i].positionBEnd+1][j]->align, category,
                                          score, couldMerge, megerScore, align);
                        if( 2 == couldMerge ){
                            break;
                        } else if (couldMerge == 1) {
                            allExtendSeeds[i].positionAEnd = queryPositionAsIndex[allExtendSeeds[i].positionBEnd +1][j]->positionAEnd;
                            queryPositionAsIndex[allExtendSeeds[i].positionBEnd + 1][j]->everUsed = true;
                            allExtendSeeds[i].positionBEnd = queryPositionAsIndex[allExtendSeeds[i].positionBEnd +1][j]->positionBEnd;
                            allExtendSeeds[i].align = align;
                            allExtendSeeds[i].score = megerScore;
                            goto queryPositionAsIndex;
                        }
                    }
                }
            }
        }
    }
    queryPositionAsIndex.clear();
    std::cout << "begin to shrink allMergedSeeds " << allExtendSeeds.size() << std::endl;
    for( i=j=0; i<allExtendSeeds.size(); ++i ){
        if( !allExtendSeeds[i].everUsed && (allExtendSeeds[i].align.size())>=miniseedLength  ){
            if( i != j ){
                allExtendSeeds[j]=allExtendSeeds[i];
            }
            ++j;
        }
    }
    allExtendSeeds.resize(j);
    my_time = time(NULL); printf("%s", ctime(&my_time));

    std::cout << "before chain allMergedSeeds size " << allExtendSeeds.size() << std::endl;
    std::ofstream ofile;
    ofile.open("/Users/bs674/WSA_Alignment");

    for( i=0; i<allExtendSeeds.size(); ++i ) {
//        ofile << allExtendSeeds[i].positionAStart << " " << allExtendSeeds[i].positionBStart <<
//                  " " << allExtendSeeds[i].positionAEnd << " " << allExtendSeeds[i].positionBEnd << " " << " " << allExtendSeeds[i].score<<
//                  allExtendSeeds[i].align << std::endl;
        //outputAlignment(seqAChar, seqBChar, allExtendSeeds[i].positionAStart, allExtendSeeds[i].positionBStart, allExtendSeeds[i].align);
    }

    std::vector<MergedSeed> chain;
    longestPath_V1 ( allExtendSeeds, chain );
    allExtendSeeds.clear();
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "final chain " << chain.size() << std::endl;

    int32_t AStart = 0;
    int32_t BStart = 0;


    std::stringstream alignStream;
    std::stack<char> A;
    for( i=0; i<chain.size(); ++i ){
//        std::cout << chain[i].positionAStart << " " << chain[i].positionBStart <<
//                  " " << chain[i].positionAEnd << " " << chain[i].positionBEnd << " " << " " << chain[i].score << std::endl << chain[i].align << std::endl;
//        outputAlignment(seqAChar, seqBChar, chain[i].positionAStart, chain[i].positionBStart, chain[i].align);

        zdp(seqAChar+AStart, seqBChar+BStart, chain[i].positionAStart-AStart, chain[i].positionBStart-BStart, category+AStart, A, score);
        while (!A.empty()) {
            alignStream << A.top();
            A.pop();
        }
        alignStream << chain[i].align;
        AStart = chain[i].positionAEnd+1;
        BStart = chain[i].positionBEnd+1;
    }
    std::cout << "after chain" << std::endl;
    zdp(seqAChar+AStart, seqBChar+BStart, lengthSeqA-AStart, lengthSeqB-BStart, category+AStart, A, score);
    while (!A.empty()) {
        alignStream << A.top();
        A.pop();
    }

    std::string thisAlign = alignStream.str();
    std::string aa;
    std::string ab;
    std::string cc;
    int ai=0;
    int bi=0;
    for (  i=0; i<thisAlign.size(); ++i ){
        if( thisAlign[i] == 'M' ){
            aa += seqA[ai];
            ab += seqB[bi];
            int16_t printableScore = category[ai]/10-1;
            cc += std::to_string(printableScore);
            ++ai;
            ++bi;
        }else if( thisAlign[i] == 'D' ){
            aa += seqA[ai];
            ab += '-';
            int16_t printableScore = category[ai]/10-1;
            cc += std::to_string(printableScore);
            ++ai;
        }else if( thisAlign[i] == 'I' ){
            aa += '-';
            ab += seqB[bi];
            cc += '-';
            ++bi;
        }
    }

    std::size_t pos = 0;
    int16_t linewidth = 180;
    int32_t currentALength = 0;
    int32_t currentBLength = 0;
    while (pos < aa.length()) {
        std::string t1 = aa.substr(pos, linewidth);
        std::string t2 = ab.substr(pos, linewidth);
        for( i=0; i<t1.size(); ++i ){
            if( t1[i] != '-' ){
                ++currentALength;
            }
            if( t2[i] != '-' ){
                ++currentBLength;
            }
        }

        ofile << aa.substr(pos, linewidth) << " " << currentALength << std::endl;
        ofile << thisAlign.substr(pos, linewidth) << std::endl;
        ofile << ab.substr(pos, linewidth) << " " << currentBLength << std::endl;
        ofile << cc.substr(pos, linewidth) << std::endl << std::endl << std::endl << std::endl;
        pos += linewidth;
    }

//    std::cout << aa << std::endl << thisAlign << std::endl << ab << std::endl << cc << std::endl;
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << std::endl;
}

*/