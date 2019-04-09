//
// Created by Baoxing Song on 2019-04-02.
//

#include <cstdint>
#include "seed.h"

/*
struct InitialSeed {
    int32_t positionAStart;
    int32_t positionBStart;
    int8_t score;
};

struct ExtendSeed {
    int32_t length;
    int32_t positionAStart;
    int32_t positionBStart;
    int64_t score;
};

struct MergedSeed {
    int32_t positionAStart;
    int32_t positionBStart;
    int32_t positionAEnd;
    int32_t positionBEnd;
    std::string align; // if align is "", it means every should be M, the base align is indicated using the cigar letters M D I
    int64_t score;
};


//O2
void creatSeed(const char * a, const char * b, const int8_t & seedLength, int8_t & score){ // never use initial seed length longer than 126
    score=0;
    for ( int8_t i=0; i < seedLength; ++i ) {
        if ( *(a+i) == *(b+i) ){
            ++score;
        }else{
            --score;
        }
    }
}

//computational complexity should be o1
void forwardExtendSeed(const char * a, const char * b, const int32_t & sequenceLengthA, const int32_t & sequenceLengthB, const int8_t & seedLength, const int8_t & seedScore, int32_t & extendedSeedLength, int32_t & extendSeedScore){
    extendedSeedLength = seedLength;
    extendSeedScore = seedScore;
    int currentScore = seedScore;
    for ( int32_t i=seedLength; i < sequenceLengthA && i < sequenceLengthB; ++i ) {
        if ( *(a+i) == *(b+i) ){
            ++currentScore;
        }else{
            return;
            --currentScore;
        }
        if ( currentScore > extendSeedScore ) {
            extendSeedScore = currentScore;
            extendedSeedLength = i + 1;
        }else if( 0 >= currentScore ){
            return;
        }
    }
}

// we are not using it
void reverseExtendSeed(const char * a, const char * b, const int32_t & positionA, const int32_t & positionB,
        const int32_t & seedLength, const int32_t & seedScore, int32_t & reverseExtendedLength, int32_t & reverseExtendSeedScore){
    reverseExtendedLength = 0;
    reverseExtendSeedScore = seedScore;
    int miniPosition = positionA<positionB ? positionA : positionB;
    if( miniPosition>0 ){
        int currentScore = seedScore;
        for ( int32_t i=1; i <=miniPosition; ++i ) {
            if ( *(a-i) == *(b-i) ){
                ++currentScore;
            }else{
                return;
                --currentScore;
            }
            if ( currentScore > reverseExtendSeedScore ) {
                reverseExtendSeedScore = currentScore;
                reverseExtendedLength = i;
            }else if( 0 >= currentScore ){
                return;
            }
        }
    }
}


//here I assume position1A < position2A and  position2B < position2B
// A means the reference seeds, B means the query seeds
// 1 means the upstream seed and 2 meads the downstream seed
// couldMeger 1 means megered
//            0 those two seeds overlapped, should check next one
//            2 the second seed is too far from the first one, should stop searching
//            3 other reasons, should check next one
//  computational complexity should be O1
void megerSeed ( const int & seedScore1, const int & seedScore2, const int32_t & position1AEnd,
        const int32_t & position2AStart, const int32_t & position1BEnd, const int32_t & position2BStart,
        const int32_t & seedLength1, const int32_t & seedLength2, const int16_t * category, Score & score,
        int8_t & couldMerge, int32_t & megerScore, std::string & align ) {
    align="";
    couldMerge = 0;
    if( (position1AEnd+1 > position2AStart) ){
        return;
    }else if ( (position1AEnd+1 < position2AStart) ) {
        if( (position1BEnd+1) < position2BStart  ){
            couldMerge=2;
            return;
        }else if ( (position1BEnd+1) > position2BStart ){
            return;
        }else if ( (position1BEnd+1) == position2BStart ) {
            int32_t gapPenalty; // here using the smaller indel penalty
            if( score.getOpenPenalty(*(category+position1AEnd+1)) > score.getOpenPenalty(*(category+position2AStart-1)) ){
                gapPenalty=score.getOpenPenalty(*(category+position2AStart-1));
                for( int i=0; position2AStart-1-i > position1AEnd; ++i){
                    gapPenalty += score.getExtendPenalty(*(category+position2AStart-1-i));
                }
            }else{
                gapPenalty=score.getOpenPenalty(*(category+position1AEnd+1));
                for( int i=0; position1AEnd+1+i< position2AStart; ++i){
                    gapPenalty += score.getExtendPenalty(*(category+position1AEnd+1+i));
                }
            }
            int32_t thisScore = seedScore1 + seedScore2 - gapPenalty;
            //int32_t thisScore = seedScore1 + seedScore2 - position2AStart + position1AEnd;
            if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
                megerScore = thisScore;
                couldMerge = 1;
            }else{
                return;
            }
            align = std::string(seedLength1, 'M')+std::string((position2AStart-(position1AEnd+1)), 'D')+std::string(seedLength2, 'M');
        }
    }else if ( (position1AEnd+1 == position2AStart)  ){
        if( (position1BEnd+1) > position2BStart  ) {
            return; // check next
        }else if( (position1BEnd+1) == position2BStart  ){
            int32_t thisScore = seedScore1 + seedScore2;
            if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
                megerScore = thisScore;
                couldMerge = 1;
            }else{
                return;
            }
            align = std::string(seedLength1+seedLength2, 'M');
            std::cerr << "two seeds could be linked together directly, this should never happen." << std::endl;

        }else if ((position1BEnd+1) < position2BStart)  {
            int16_t openGapPenalty, extendGapPenalty; // here using the smaller indel penalty
            if( score.getOpenPenalty(*(category+position1AEnd)) > score.getOpenPenalty(*(category+position2AStart)) ){
                openGapPenalty=score.getOpenPenalty(*(category+position2AStart));
                extendGapPenalty=score.getExtendPenalty(*(category+position2AStart));
            }else{
                openGapPenalty=score.getOpenPenalty(*(category+position1AEnd));
                extendGapPenalty=score.getExtendPenalty(*(category+position1AEnd));
            }

            int32_t thisScore = seedScore1 + seedScore2 - openGapPenalty - (position2BStart-(position1BEnd+1))*extendGapPenalty; //extend gap penalty  // assume we have a seed 010-010 the first seed score is 2 and the second seed score is 2 and the merged seed score is 2 (0 means match, 1 means mis-match and - means indel)
            //int32_t thisScore = seedScore1 + seedScore2 - position2BStart + position1BEnd;
            if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
                megerScore = thisScore;
                couldMerge = 1;
            }else{
                return;
            }
            align = std::string(seedLength1, 'M') + std::string((position2BStart-(position1BEnd+1)), 'I') + std::string(seedLength2, 'M');
        }
    }else{
        std::cerr << "megerSeed, this line should never run" << std::endl;
    }
}

//o1 but the number of seeds is o2 so acctually yhe time costing is o2
void megerSeed ( const int & seedScore1, const int & seedScore2, const int32_t & position1AEnd,
                 const int32_t & position2AStart, const int32_t & position1BEnd, const int32_t & position2BStart,
                 const std::string & align1, const std::string & align2, const int16_t * category,
                 Score & score, int8_t & couldMerge, int32_t & newScore,
                 std::string & align ) {
    couldMerge = 0;
    if( position1AEnd >= position2AStart ){
        return;
    }else if ( (position1AEnd+1 < position2AStart) && (position1BEnd+1) < position2BStart ) {
        couldMerge=2;
        return;
    }else if ( (position1AEnd+1 == position2AStart) && (position1BEnd+1) < position2BStart ){

        int16_t openGapPenalty, extendGapPenalty; // here using the smaller indel penalty
        if( score.getOpenPenalty(*(category+position1AEnd)) > score.getOpenPenalty(*(category+position2AStart)) ){
            openGapPenalty=score.getOpenPenalty(*(category+position2AStart));
            extendGapPenalty=score.getExtendPenalty(*(category+position2AStart));
        }else{
            openGapPenalty=score.getOpenPenalty(*(category+position1AEnd));
            extendGapPenalty=score.getExtendPenalty(*(category+position1AEnd));
        }

        int32_t thisScore = seedScore1 + seedScore2 - openGapPenalty - (position2BStart-(position1BEnd+1))*extendGapPenalty;  // extend gap penalty assume we have a seed 010-010 the first seed score is 2 and the second seed score is 2 and the merged seed score is 2 (0 means match, 1 means mis-match and - means indel)
        //int32_t thisScore = seedScore1 + seedScore2 - position2BStart + position1BEnd;
        if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
            newScore = thisScore;
            couldMerge = 1;
        }
        align = align1 + std::string((position2BStart-(position1BEnd+1)), 'I') + align2;
    }else if( (position1AEnd+1 < position2AStart) && (position1BEnd+1) == position2BStart ){
        int32_t gapPenalty; // here using the smaller indel penalty
        if( score.getOpenPenalty(*(category+position1AEnd+1)) > score.getOpenPenalty(*(category+position2AStart-1)) ){
            gapPenalty=score.getOpenPenalty(*(category+position2AStart-1));
            for( int i=0; position2AStart-1-i > position1AEnd; ++i){
                gapPenalty += score.getExtendPenalty(*(category+position2AStart-1-i));
            }
        }else{
            gapPenalty=score.getOpenPenalty(*(category+position1AEnd+1));
            for( int i=0; position1AEnd+1+i< position2AStart; ++i){
                gapPenalty += score.getExtendPenalty(*(category+position1AEnd+1+i));
            }
        }
        int32_t thisScore = seedScore1 + seedScore2 - gapPenalty;
        //int32_t thisScore = seedScore1 + seedScore2 - position2AStart + position1AEnd;
        if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
            newScore = thisScore;
            couldMerge = 1;
        }
        align = align1+ std::string((position2AStart-(position1AEnd+1)), 'D') +align2;
    }else if ( (position1AEnd+1 == position2AStart) && (position1BEnd+1) == position2BStart ){
        int32_t thisScore = seedScore1 + seedScore2;
        if ( thisScore >= seedScore1 && thisScore >= seedScore2 ){
            newScore = thisScore;
            couldMerge = 1;
        }
        align = align1+align2;
        std::cout << "two seeds could be linked together directly, this should never happen." << std::endl;
    }else{
        couldMerge=3;
    }
}


// the computational complexity of this step should be o1
// only apply this function on seed without indel
void removeExtendSeedsDuplications(std::vector<ExtendSeed> & allExtendSeeds){
    std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](ExtendSeed a, ExtendSeed b) {
        return a.positionAStart < b.positionAStart;
    });
    int i, j;
    for ( i=0, j=0; i< allExtendSeeds.size(); ++i ){
        if( allExtendSeeds[j].positionAStart <= allExtendSeeds[i].positionAStart && ( (allExtendSeeds[j].positionAStart+allExtendSeeds[j].length) >= (allExtendSeeds[i].positionAStart+allExtendSeeds[i].length) )
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
void longestPath (std::vector<MergedSeed> & allMergedSeeds, std::vector<MergedSeed> & chain){
    std::sort(allMergedSeeds.begin(), allMergedSeeds.end(), [](MergedSeed a, MergedSeed b) {
        return a.positionAStart < b.positionAStart;
    });
    double maxSore = 0;
    int bestEnd = 0;
    double scoreArray [allMergedSeeds.size()]; // arrays of scores
    int prev [allMergedSeeds.size()];  // index of previous node in longest path
    scoreArray[0] = allMergedSeeds[0].score;
    prev[0] = -1;

    if (scoreArray[0] > maxSore){
        bestEnd = 0;
        maxSore = scoreArray[0];
    }

    for (int idx = 1; idx < allMergedSeeds.size(); ++idx) {
        scoreArray[idx] = allMergedSeeds[idx].score;
        prev[idx] = -1;
        for (int jdx = idx - 1; jdx >= 0; --jdx) {// checking all previous nodes
            // Because we swapped asm/query start position so that inversions were all increasing,
            // we should always be on the diagonal.  If not, then we filter it.
            // This gets rid of the noise, while preserving the inversions on
            // the diagonal
            // Are only looking at positions previous to our current "idx" position
            if ( (scoreArray[jdx] + allMergedSeeds[idx].score) > scoreArray[idx] &&
                allMergedSeeds[jdx].positionAEnd < allMergedSeeds[idx].positionAStart &&
                allMergedSeeds[jdx].positionBEnd < allMergedSeeds[idx].positionBStart ){
                scoreArray[idx] = scoreArray[jdx] + allMergedSeeds[idx].score;
                prev[idx] = jdx;
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

    int32_t lengthSeqA = seqA.size();
    int32_t lengthSeqB = seqB.size();
    //std::cout << "lengthSeqA: " << lengthSeqA << std::endl;
    const char * seqAChar = seqA.c_str();
    const char * seqBChar = seqB.c_str();
    std::vector<InitialSeed> allInitialSeeds;
    int8_t initialSeedsScore;
    int32_t i, j, k, l;
    for( i=0; i<=(lengthSeqA-initialSeedLength); ++i ){
        for( j=0; j<=(lengthSeqB-initialSeedLength); ++j ){
            creatSeed(seqAChar+i, seqBChar+j, initialSeedLength, initialSeedsScore);
            if( initialSeedsScore > initialSeedsScoreThreadsHold ){
                InitialSeed seed = { i, j, initialSeedsScore };
                allInitialSeeds.push_back(seed);
            }
        }
    }
    std::cout << "Initial Seeds" << std::endl;
    for( i=0; i<allInitialSeeds.size(); ++i ){
        //std::cout << allInitialSeeds[i].positionAStart << " " << allInitialSeeds[i].positionBStart << " " << allInitialSeeds[i].score << std::endl;
    }

    // seed extend begin
    int32_t extendedSeedLength;
    int32_t extendSeedScore;
    std::vector<ExtendSeed> allExtendSeeds;
    for( i=0; i<allInitialSeeds.size(); ++i ){
        InitialSeed seed = allInitialSeeds[i];
        forwardExtendSeed( seqAChar+seed.positionAStart, seqBChar+seed.positionBStart, lengthSeqA-seed.positionAStart,
                lengthSeqB-seed.positionBStart, initialSeedLength, seed.score, extendedSeedLength, extendSeedScore);
        ExtendSeed extendSeed={extendedSeedLength, seed.positionAStart, seed.positionBStart, extendSeedScore};
        allExtendSeeds.push_back(extendSeed);
    }
    allInitialSeeds.clear(); // RAM saving
    std::cout << "forward extend of seeds" << std::endl;
//    for( i=0; i<allExtendSeeds0.size(); ++i ){
//        std::cout << allExtendSeeds0[i].positionAStart << " " << allExtendSeeds0[i].positionBStart << " " << allExtendSeeds0[i].score << " " << allExtendSeeds0[i].length << std::endl;
//    }
    removeExtendSeedsDuplications(allExtendSeeds); // remove duplication records
    std::cout << "forward extend of seeds remove duplication" << std::endl;

    for( i=0; i<allExtendSeeds.size(); ++i ){
        int64_t newScore=0;
        for( j=0; j< allExtendSeeds[i].length; ++j){
            newScore += score.getScore( *(category+j+allExtendSeeds[i].positionAStart), *(seqAChar+j+allExtendSeeds[i].positionAStart), *(seqBChar+j+allExtendSeeds[i].positionBStart) );
        }
        allExtendSeeds[i].score=newScore;
    }
    std::cout << "update seed score according to their category " << std::endl;
    // update seed score according to their category end
    time_t my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << std::endl;
    // the following code could be dangerous, should be really careful
    // since there are INDELs, so should not remove duplication by comparing the start and end position without consider the alignment
    std::vector<MergedSeed> allMergedSeeds(allExtendSeeds.size());
    int8_t couldMerge=2;
    std::string align;
    for( i=0; i<allExtendSeeds.size(); ++i ){
        std::string align (allExtendSeeds[i].length, 'M');
        MergedSeed mergedSeed={allExtendSeeds[i].positionAStart, allExtendSeeds[i].positionBStart,allExtendSeeds[i].positionAStart+allExtendSeeds[i].length-1, allExtendSeeds[i].positionBStart+allExtendSeeds[i].length-1, align, allExtendSeeds[i].score};
        allMergedSeeds[i]=(mergedSeed);
    }
    allExtendSeeds.clear(); //RAM saving
    std::sort(allMergedSeeds.begin(), allMergedSeeds.end(), [](MergedSeed a, MergedSeed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });

    std::cout << "first around megered seeds" << std::endl;
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << std::endl;
//    for( i=0; i<allMergedSeeds0.size(); ++i ){
//        std::cout << allMergedSeeds0[i].positionAStart << " " << allMergedSeeds0[i].positionBStart <<
//                  " " << allMergedSeeds0[i].positionAEnd << " " << allMergedSeeds0[i].positionBEnd << " " << allMergedSeeds0[i].align << " " << allMergedSeeds0[i].score << std::endl;
//    }
    int32_t megerScore;

    size_t allMergedSeedsSize = allMergedSeeds.size();
    int32_t megerSeedStepWidth=0;
    if(allMergedSeedsSize>100){
        std::vector<int32_t> endPoints(100);
        for( i=0; i<100; ++i ){
                for (j = i + 1; j < allMergedSeedsSize; ++j) {
                    megerSeed(allMergedSeeds[i].score, allMergedSeeds[j].score, allMergedSeeds[i].positionAEnd,
                              allMergedSeeds[j].positionAStart, allMergedSeeds[i].positionBEnd,
                              allMergedSeeds[j].positionBStart,
                              allMergedSeeds[i].align, allMergedSeeds[j].align, category,
                              score, couldMerge, megerScore, align);
                    if (couldMerge == 1) { // very good
                        MergedSeed mergedSeed;
                        mergedSeed.positionAStart = allMergedSeeds[i].positionAStart;
                        mergedSeed.positionBStart = allMergedSeeds[i].positionBStart;
                        mergedSeed.positionAEnd = allMergedSeeds[j].positionAEnd;
                        mergedSeed.positionBEnd = allMergedSeeds[j].positionBEnd;
                        mergedSeed.align = align;
                        mergedSeed.score = megerScore;
                        allMergedSeeds[i] = mergedSeed;
                        --i;
                        break;
                    } else if (couldMerge == 2) {
                        break;
                    }
                }
            if( i<100 ){
                endPoints[i]=j;
            }
        }
        std::sort(endPoints.begin(), endPoints.end(), [](int32_t a, int32_t b) {
            return a < b;
        });
        megerSeedStepWidth=endPoints[25];
    }else{
        std::cerr << "the sequence is too short or the seed length is to long to analysis using this approach" << std::endl;
        exit(1);
    }
    std::set<int32_t> usedIndex;
    int32_t thisMegerSeedStepWidth;
    for( i=0; i<allMergedSeedsSize; ++i ){
        std::cout << "i " << i << " of " << allMergedSeedsSize << std::endl;
        if( usedIndex.find(i) == usedIndex.end() ) {
            j = i + 1;
            // for speeding up begin
            thisMegerSeedStepWidth=megerSeedStepWidth;
            k=j+thisMegerSeedStepWidth;
            while( (allMergedSeeds[i].positionAEnd) >= allMergedSeeds[k].positionAStart && thisMegerSeedStepWidth>5){
                j=k;
                k+=thisMegerSeedStepWidth;
                if( k>=allMergedSeedsSize ){
                    break;
                }
                while( (allMergedSeeds[i].positionAEnd) < allMergedSeeds[k].positionAStart && (allMergedSeeds[i].positionBEnd) < allMergedSeeds[k].positionBStart && thisMegerSeedStepWidth>1 ){
                    thisMegerSeedStepWidth=thisMegerSeedStepWidth/2;
                    k=j+thisMegerSeedStepWidth;
                }
            }
            // for speeding up end
            for (; j < allMergedSeedsSize; ++j) {
                if( usedIndex.find(j)==usedIndex.end() ) {
                    megerSeed(allMergedSeeds[i].score, allMergedSeeds[j].score, allMergedSeeds[i].positionAEnd,
                              allMergedSeeds[j].positionAStart, allMergedSeeds[i].positionBEnd,
                              allMergedSeeds[j].positionBStart,
                              allMergedSeeds[i].align, allMergedSeeds[j].align, category,
                              score, couldMerge, megerScore, align);
                    if (couldMerge == 1) { // very good
                        MergedSeed mergedSeed;
                        mergedSeed.positionAStart = allMergedSeeds[i].positionAStart;
                        mergedSeed.positionBStart = allMergedSeeds[i].positionBStart;
                        mergedSeed.positionAEnd = allMergedSeeds[j].positionAEnd;
                        mergedSeed.positionBEnd = allMergedSeeds[j].positionBEnd;
                        mergedSeed.align = align;
                        mergedSeed.score = megerScore;
                        allMergedSeeds[i] = mergedSeed;
                        --i;
                        usedIndex.insert(j);
                        break;
                    } else if (couldMerge == 2) {
                        break;
                    }
                }
            }
        }
    }

    for( i=j=0; i<allMergedSeeds.size(); ++i ){
        if( usedIndex.find(i)==usedIndex.end() && (allMergedSeeds[i].positionAEnd-allMergedSeeds[i].positionAStart)>=miniseedLength  ){ // should be kept
            if( i != j ){
                allMergedSeeds[j]=allMergedSeeds[i];
            }
            ++j;
        }
    }
    allMergedSeeds.resize(j);
    std::cout << "final merged seeds" << std::endl;
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << std::endl;

    for( i=0; i<allMergedSeeds.size(); ++i ){
        std::cout << allMergedSeeds[i].positionAStart << " " << allMergedSeeds[i].positionBStart <<
            " " << allMergedSeeds[i].positionAEnd << " " << allMergedSeeds[i].positionBEnd << " align " << allMergedSeeds[i].align << " " << allMergedSeeds[i].score << std::endl;
    }

    std::cout << "before chain allMergedSeeds size " << allMergedSeeds.size() << std::endl;
    std::vector<MergedSeed> chain;
    longestPath ( allMergedSeeds, chain );
    allMergedSeeds.clear();
    std::cout << "final chain" << std::endl;
    for( i=0; i<chain.size(); ++i ){
        std::cout << chain[i].positionAStart << " " << chain[i].positionBStart <<
                  " " << chain[i].positionAEnd << " " << chain[i].positionBEnd << " " << chain[i].align << " " << chain[i].score << std::endl;
        outputAlignment(seqAChar, seqBChar, chain[i].positionAStart, chain[i].positionBStart, chain[i].align);
    }
}
*/
s