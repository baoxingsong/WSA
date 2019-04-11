/*
 * =====================================================================================
 *
 *       Filename:  SuffixTree.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  08/04/2019 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, songbaoxing168@163.com
 *
 * =====================================================================================
 */

/*************************************************************************

Each non-root node has one and only one parent node.
Each non-leaf node has <=4 son nodes.
Each node was created only when the sequence has been observed.
Each leaf node could tell us the position if this branch on the reference sequence
Each node is a char from the reference genome sequence (here we set it as not string)
Since we have the a ROOT so the level of the tree is k+1 (k is the k-mer size)
The level of leaves is always k+!

 ************************************************************************/


#ifndef WSA_SUFFIXTREE_H
#define WSA_SUFFIXTREE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

struct NODE {
    /* A pointer to that node's father */
    int32_t father;
    /* the DNA sequence charactor of this NODE */
    char c;
    /* which level of this seed is located on the tree*/
    int8_t level;
    /* A linked list of sons of that node */
    // the index in the node vector
    std::vector<int32_t> sons;
    // this is only for leaves
    // which is the position of this branch on the reference sequence
    std::vector<int32_t> positons;
};

void CreateTree(const char * seq, const int32_t & seq_length, const int8_t & k, std::vector<NODE> & allNodes);

void findSubString( std::vector<NODE> & allNodes, const char * qSeq, const int32_t & qSeqLength, const int8_t & k, const int8_t & maximumMismatch,
                    std::vector<int32_t > & databasePositions, std::vector<int32_t > & queryPositions, std::vector<int8_t> & numMisMacth );
void findSubString( std::vector<NODE> & allNodes, const char * qSeq, const int8_t & k, const int8_t & maximumMismatch,
                    std::vector<int32_t > & databasePositions, std::vector<int32_t > & queryPositions, std::vector<int8_t> & numMisMacth );

#endif //WSA_SUFFIXTREE_H
