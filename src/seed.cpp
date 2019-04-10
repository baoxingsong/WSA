//
// Created by Baoxing Song on 2019-04-02.
//

#include <cstdint>
#include "seed.h"

struct Seed {
    int32_t positionAStart;
    int32_t positionBStart;
    int64_t score;

    int32_t positionAEnd;
    int32_t positionBEnd;
    std::string align; // if align is "", it means every should be M, the base align is indicated using the cigar letters M D I
    bool everUsed;
};

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
        int32_t thisScore = seedScore1 + seedScore2 - openGapPenalty/*open gap penalty*/ - (position2BStart-(position1BEnd+1))*extendGapPenalty/*extend gap penalty*/;  // assume we have a seed 010-010 the first seed score is 2 and the second seed score is 2 and the merged seed score is 2 (0 means match, 1 means mis-match and - means indel)
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
// this function does not work so well as expected

void removeExtendSeedsDuplications_V_new(std::vector<Seed> & allExtendSeeds){

    // this part was implemented basing the assumption that the seeds, should be merged, share the same end position
    std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](Seed a, Seed b) {
        if(a.positionAEnd == b.positionAEnd){
            if( a.positionBEnd == b.positionBEnd){
                return a.positionAStart < b.positionAStart;
            }else{
                return a.positionBEnd < b.positionBEnd;
            }
        }else{
            return a.positionAEnd < b.positionAEnd;
        }
    });
//    std::cout << "line 148" << std::endl;
    int i, j;
    for ( i=0, j=0; i< allExtendSeeds.size(); ++i ){
//        if( 0 == i% 1000 ){
//            std::cout << "line 152 i " << i << std::endl;
//        }
        if(  ( (allExtendSeeds[j].positionAEnd) == (allExtendSeeds[i].positionAEnd) ) && ( (allExtendSeeds[j].positionBEnd) == (allExtendSeeds[i].positionBEnd ))
             && allExtendSeeds[j].positionAStart <= allExtendSeeds[i].positionAStart && allExtendSeeds[j].positionBStart <= allExtendSeeds[i].positionBStart  ){

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
void longestPath (std::vector<Seed> & allMergedSeeds, std::vector<Seed> & chain){
    std::sort(allMergedSeeds.begin(), allMergedSeeds.end(), [](Seed a, Seed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    std::cout << "line 161 begin to chain" << std::endl;
    int64_t maxSore = 0;
    int32_t bestEnd = -1;
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

void longestPath (std::vector<Seed> & allMergedSeeds, int32_t & startIndex, const int32_t & referenceStart,
        const int32_t & referenceEnd, const int32_t & queryStart, const int32_t & queryEnd, std::vector<Seed> & chain){
    int64_t maxSore = 0;
    int32_t bestEnd = -1;
    std::vector<int64_t> scoreArray (allMergedSeeds.size()); // arrays of scores
    std::vector<int32_t> prev (allMergedSeeds.size());  // index of previous node in longest path
    scoreArray[startIndex] = allMergedSeeds[startIndex].score;
    prev[startIndex] = -1;
    if( allMergedSeeds[startIndex].positionAStart>referenceStart && allMergedSeeds[startIndex].positionAEnd < referenceEnd &&
        allMergedSeeds[startIndex].positionBStart>queryStart && allMergedSeeds[startIndex].positionBEnd < queryEnd ) {
        if (scoreArray[startIndex] > maxSore) {
            bestEnd = startIndex;
            maxSore = scoreArray[startIndex];
        }
    }
    for (int32_t idx = startIndex+1; idx < allMergedSeeds.size(); ++idx) {
        scoreArray[idx] = allMergedSeeds[idx].score;
        prev[idx] = -1;
        if( allMergedSeeds[idx].positionAStart >= referenceEnd ){
            break;
        }
        if( allMergedSeeds[idx].positionAStart > referenceStart && allMergedSeeds[idx].positionAEnd < referenceEnd &&
            allMergedSeeds[idx].positionBStart > queryStart && allMergedSeeds[idx].positionBEnd < queryEnd ){
            for (int32_t jdx = idx - 1; jdx >= startIndex; --jdx) {// checking all previous nodes
                // Are only looking at positions previous to our current "idx" position
                if( allMergedSeeds[jdx].positionAStart <= referenceStart ){
                    break;
                }
                if (allMergedSeeds[jdx].positionAEnd < allMergedSeeds[idx].positionAStart &&
                    allMergedSeeds[jdx].positionBEnd < allMergedSeeds[idx].positionBStart &&
                    allMergedSeeds[jdx].positionAStart>referenceStart && allMergedSeeds[jdx].positionAEnd < referenceEnd &&
                    allMergedSeeds[jdx].positionBStart > queryStart && allMergedSeeds[jdx].positionBEnd < queryEnd ) {
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
    }
    if( bestEnd>=0 ){
        int idx=bestEnd; // bestEnd is where to stop the longest path
        chain.push_back(allMergedSeeds[idx]);
        int jdx = prev[idx]; // prev[] is index on the longest path
        while( jdx>=0 ){
            chain.push_back(allMergedSeeds[jdx]);
            jdx=prev[jdx];
        }// Reversing the order
        startIndex=bestEnd+1;
    }
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



void seq2seed ( const char * seqAChar, const char * seqBChar, const int32_t & lengthSeqA, const int32_t & lengthSeqB, const int16_t * category,
                const int8_t & initialSeedLength, const int8_t & initialSeedsScoreThreadsHold,
                const std::string & outputFile, int16_t miniseedLength, Score & score,
                const int64_t & numberOfSeedsForChain){
    miniseedLength--;

    int8_t initialSeedsScore;
    int32_t i, j;
    int32_t extendedSeedLength;
    int32_t extendSeedScore;
    std::vector<Seed> allExtendSeeds;

    std::vector<NODE> allNodes;
    CreateTree(seqAChar, lengthSeqA, initialSeedLength, allNodes);
    int8_t maximumMismatch=initialSeedLength-initialSeedsScoreThreadsHold;
    std::vector<int32_t > databasePositions;
    std::vector<int32_t > queryPositions;
    std::vector<int8_t > numMisMacth;
    std::cout << "line 323" << std::endl;
    for(  i=0; i<=(lengthSeqB-initialSeedLength); ++i) {
        databasePositions.clear();
        queryPositions.clear();
        numMisMacth.clear();
        findSubString(allNodes, seqBChar+i, initialSeedLength, maximumMismatch, databasePositions, queryPositions, numMisMacth);
        for ( j = 0; j < databasePositions.size(); ++j) {
            initialSeedsScore = initialSeedLength-numMisMacth[j];
            forwardExtendSeed( seqAChar+databasePositions[j], seqBChar+i, lengthSeqA-databasePositions[j],
                               lengthSeqB-i, initialSeedLength, initialSeedsScore, extendedSeedLength, extendSeedScore);
            Seed extendSeed={databasePositions[j], i, extendSeedScore, databasePositions[j]+extendedSeedLength-1,
                             i+extendedSeedLength-1, std::string(extendedSeedLength, 'M'), false};
            allExtendSeeds.push_back(extendSeed);
        }
        if( i>0 && allExtendSeeds.size()>1000 && 0 == i%50000 ){
//            std::cout << "line 374 i " << i << std::endl;
            removeExtendSeedsDuplications_V_new(allExtendSeeds);
//            std::cout << "line 376 i " << i << std::endl;
        }
    }

    std::cout << "Initial Seeds " << allExtendSeeds.size() << std::endl;
    removeExtendSeedsDuplications_V_new(allExtendSeeds); // remove duplication records
    std::cout << " seeds remove duplication " << allExtendSeeds.size() << std::endl;

    // update seed score according to their category
    int64_t newScore;
    for( i=0; i<allExtendSeeds.size(); ++i ){
        newScore=0;
        for( j=allExtendSeeds[i].positionAStart; j< allExtendSeeds[i].positionAEnd; ++j){
            newScore += score.getScore( *(category+j), *(seqAChar+j), *(seqBChar+j-allExtendSeeds[i].positionAStart+allExtendSeeds[i].positionBStart) ); // the category share the same coordinate with reference sequence
        }
        allExtendSeeds[i].score=newScore;
    }
    time_t my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update seed score according to their category " << allExtendSeeds.size() << std::endl;

    // the following code could be dangerous, should be really careful
    // since there are INDELs, so should not remove duplication by comparing the start and end position without consider the alignment
    //std::vector<MergedSeed> allMergedSeeds(allExtendSeeds.size());

    std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](Seed a, Seed b) {
        if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
    });
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "extend transformed to megered seeds " << allExtendSeeds.size() << std::endl;

    // prepare special data structure for speeding up purpose
    std::vector<std::vector<Seed *>> referencePositionAsIndex(lengthSeqA+1); // plus for end based entry
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
    ++j;
    allExtendSeeds.resize(j);
    std::cout << "referencePositionAsIndex shrink allMergedSeeds " << allExtendSeeds.size() << std::endl;
    allMergedSeedsSize=allExtendSeeds.size();
    std::vector<std::vector<Seed *>> queryPositionAsIndex(lengthSeqB+1);
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
    ++j;
    allExtendSeeds.resize(j);
    my_time = time(NULL); printf("%s", ctime(&my_time));

    std::cout << "before chain allMergedSeeds size " << allExtendSeeds.size() << std::endl;
    std::ofstream ofile;
    ofile.open(outputFile);

    std::vector<Seed> chain;
    if( numberOfSeedsForChain>0 && allExtendSeeds.size() > numberOfSeedsForChain ){
        std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](Seed a, Seed b) {
            return  a.score > b.score;
        });
        std::vector<Seed> subMergedSeedVector(numberOfSeedsForChain);
        for( i=0; i< numberOfSeedsForChain; ++i  ){
            subMergedSeedVector[i]=allExtendSeeds[i];
        }
        longestPath ( subMergedSeedVector, chain );
        std::sort(allExtendSeeds.begin(), allExtendSeeds.end(), [](Seed a, Seed b) {
            if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
        });

        int32_t referenceStart=-1;
        int32_t queryStart=-1;
        int32_t startIndex=0;
        int32_t referenceEnd;
        int32_t queryEnd;

        size_t chainSize = chain.size();
        for( i=0; i<chainSize; ++i  ){
            //std::cout << " line 462 i " << i << " of chainSize " << chainSize << std::endl;
            referenceEnd=chain[i].positionAStart;
            queryEnd=chain[i].positionBStart;
            longestPath (allExtendSeeds, startIndex, referenceStart, referenceEnd, queryStart, queryEnd, chain);
            referenceStart=chain[i].positionAEnd;
            queryStart=chain[i].positionBEnd;
        }
        if( startIndex < allExtendSeeds.size() ){
            referenceEnd=lengthSeqA;
            queryEnd=lengthSeqB;
            longestPath (allExtendSeeds, startIndex, referenceStart, referenceEnd, queryStart, queryEnd, chain);
        }
        std::sort(chain.begin(), chain.end(), [](Seed a, Seed b) {
            if(a.positionAStart == b.positionAStart){ return a.positionBStart < b.positionBStart; }else{return a.positionAStart < b.positionAStart;}
        });
    }else{
        longestPath ( allExtendSeeds, chain );
    }

    allExtendSeeds.clear();
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "final chain " << chain.size() << std::endl;

    int32_t AStart = 0;
    int32_t BStart = 0;

    std::stringstream alignStream;
    std::stack<char> A;
    for( i=0; i<chain.size(); ++i ){
        zdp(seqAChar+AStart, seqBChar+BStart, chain[i].positionAStart-AStart, chain[i].positionBStart-BStart, category+AStart, A, score);
        while (!A.empty()) {
            alignStream << A.top();
            //std::cout << A.top();
            A.pop();
        }
        //std::cout << std::endl;
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
            aa += seqAChar[ai];
            ab += seqBChar[bi];
            int16_t printableScore = category[ai]/10-1;
            cc += std::to_string(printableScore);
            ++ai;
            ++bi;
        }else if( thisAlign[i] == 'D' ){
            aa += seqAChar[ai];
            ab += '-';
            int16_t printableScore = category[ai]/10-1;
            cc += std::to_string(printableScore);
            ++ai;
        }else if( thisAlign[i] == 'I' ){
            aa += '-';
            ab += seqBChar[bi];
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

    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << std::endl;

}
