//
// Created by Baoxing Song on 2019-04-08.
//

#include "SuffixTree.h"

// this for CreateTree
int32_t find_son( int32_t & fatherNode, const char & character, std::vector<NODE> & allNodes){
    if( allNodes[fatherNode].sons.size() == 0 ){
        return 0;
    }
    for( int32_t node : allNodes[fatherNode].sons ){
        if ( allNodes[node].c == character ) { // the current node has a son node with the specific cha
            return node;
        }
    }
    return 0;
}

void CreateTree(const char * seq, const int32_t & seq_length, const int8_t & k, std::vector<NODE> & allNodes) {
    NODE root = {-1, '0', 0};
    // the root tree is actually an empty NODE, The porpose is tell us where to start search the tree
    root.sons.resize(0);
    allNodes.push_back(root); // the root NODE is the 0 NODE
    int8_t j;
    for( int32_t i=0; i<(seq_length-k); ++i ){
        int32_t father = 0;
        for( j=0; j<k; ++j ){
            int32_t findSon = find_son(father, seq[i+j], allNodes);
            if( allNodes[findSon].level == 0 ){
                // if could not find it in the current tree, the create a new NODE and add it to the tree
                NODE thisNode = {father, seq[i+j], j+1};
                thisNode.sons.resize(0);
                allNodes.push_back(thisNode);
                allNodes[father].sons.push_back(allNodes.size()-1);
                father = allNodes.size()-1;
            }else{
                father = findSon;
            }
        }
        allNodes[father].positons.push_back(i);
    }
}

// here I am using a recursion method to search the tree
// the nodePointer is the index of the current NODE in allNodes
// thisNumberOfMismatch is the current number of mismatch, it could not be a pointer of reference
// maximumMismatch is a parameter, which is the maximum number of mismatch allowed for each seed
// queryPosition is the query position of the current k-mer sequence from the query genome sequence
// queryPositions is the position of seeds from the query genome sequence
// databasePositions is the position of seeds from the reference genome sequence
// numMisMacth is the number of mis-match of the current seed
void findSubString_p( const int32_t & nodePointer, std::vector<NODE> & allNodes, const char * qSeq, const int8_t & k, int8_t thisNumberOfMismatch,
        const int8_t & maximumMismatch, const int32_t & queryPosition, std::vector<int32_t > & databasePositions,
        std::vector<int32_t > & queryPositions, std::vector<int8_t> & numMisMacth){

    if ( allNodes[nodePointer].level == k && qSeq[k-1]==allNodes[nodePointer].c){
        // here is the leaf
        // the last char much equal
        for( int32_t p : allNodes[nodePointer].positons){
            databasePositions.push_back(p);
            numMisMacth.push_back(thisNumberOfMismatch);
            queryPositions.push_back(queryPosition);
        }
        return;
    }
    if( allNodes[nodePointer].c == qSeq[allNodes[nodePointer].level-1] ){

    }else{
        ++thisNumberOfMismatch;
    }
    if( thisNumberOfMismatch > maximumMismatch ){ // too many mismatch found on the current branch, so stop it
        return;
    }
    for( int32_t nn : allNodes[nodePointer].sons ){ // go to check next level of the tree
        findSubString_p( nn, allNodes, qSeq, k, thisNumberOfMismatch, maximumMismatch, queryPosition, databasePositions, queryPositions, numMisMacth);
    }
}

// find all the matches for a query sequence
void findSubString( std::vector<NODE> & allNodes, const char * qSeq, const int32_t & qSeqLength, const int8_t & k, const int8_t & maximumMismatch,
                    std::vector<int32_t > & databasePositions, std::vector<int32_t > & queryPositions, std::vector<int8_t> & numMisMacth ){
    for( int32_t i=0; i<=(qSeqLength-k); ++i){
        for( int32_t n : allNodes[0].sons ){
            if( allNodes[n].c == qSeq[i] ){ // the first char must equal with
                int8_t thisNumberOfMismatch = 0;
                for( int32_t nn : allNodes[n].sons ){
                    findSubString_p( nn, allNodes, qSeq+i, k, thisNumberOfMismatch, maximumMismatch, i, databasePositions, queryPositions, numMisMacth);
                }
            }
        }
    }
}

// find all the matches for the k-mer size query sequence
void findSubString( std::vector<NODE> & allNodes, const char * qSeq, const int8_t & k, const int8_t & maximumMismatch,
                    std::vector<int32_t > & databasePositions, std::vector<int32_t > & queryPositions, std::vector<int8_t> & numMisMacth ){
    int32_t i=0;
    for( int32_t n : allNodes[0].sons ){
        if( allNodes[n].c == qSeq[i] ){ // the first char must equal with
            int8_t thisNumberOfMismatch = 0;
            for( int32_t nn : allNodes[n].sons ){
                findSubString_p( nn, allNodes, qSeq+i, k, thisNumberOfMismatch, maximumMismatch, i, databasePositions, queryPositions, numMisMacth);
            }
        }
    }
}
