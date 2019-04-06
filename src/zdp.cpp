//
// Created by Baoxing Song on 2019-04-04.
//

#include "zdp.h"


/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:39
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), songbaoxing168@163.com
 *
 * =====================================================================================
 */
#include "zdp.h"
/*************************************************************************




 ************************************************************************/





void zdp(const char * seqAChar, const char * seqBChar, const int32_t & lengthA, const int32_t & lengthB, const int16_t * category, std::stack<char> & A, Score & score){
    std::cout << "lengthA " << lengthA << " lengthB " << lengthB << std::endl;
    if( lengthA==0 || 0==lengthB ){
        int32_t i = lengthB;
        int32_t j = lengthA;
        while (i > 0 || j > 0) {
            if (i == 0) {
                A.push('D');
                --j;
            } else if (j == 0) {
                A.push('I');
                --i;
            }
        }
        return;
    }

    auto ** _similarity_matrix = new int32_t*[lengthB + 1];
    for (int i = 0; i < (lengthB + 1); i++) {
        _similarity_matrix[i] = new int32_t[lengthA + 1];
    }
    for (int i = 0; i<= lengthB; i++) {
        for (int j = 0; j<= 1; j++) {
            _similarity_matrix[i][j] = 0;
        }
    }
    for (int j = 0; j<= lengthA; j++){
        for (int i = 0; i<= 1; i++) {
            _similarity_matrix[i][j] = 0;
        }
    }
    // this matrix is for set different penalty for open gap and extend gap begin
    // 0 for match, 1 for deletion, 2 for insertation
    // and the track also changed to use this matrix
    auto **_track_matrix = new VARIANTCATEGORY*[lengthB + 1];
    for (int i = 0; i < (lengthB + 1); i++) {
        _track_matrix[i] = new VARIANTCATEGORY[lengthA + 1];
    }

    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= lengthA; j++) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }
    for (int j = 0; j <= 1; j++){
        for (int i = 0; i <= lengthB; i++) {
            _track_matrix[i][j] = SNPORINSERTIONORDELETION;
        }
    }


    int l;
    int match = 0, insert = 0, del = 0;
    for (int j=1; j < lengthA + 1; ++j){
        for (l = 1; l < lengthB + 1; ++l) {
            match = _similarity_matrix[l - 1][j - 1] + score.getScore(*(category+(j-1)), *(seqAChar+j-1), *(seqBChar+l-1));

            if (l==1 || _track_matrix[l-1][j] == INSERTION || _track_matrix[l-1][j] == SNPORINSERTION
                || _track_matrix[l-1][j] == SNPORINSERTIONORDELETION || _track_matrix[l - 1][j] == INSERTIONORDELETION ) { //deletion
                insert = _similarity_matrix[l-1][j] - score.getExtendPenalty(*(category+(lengthA-1))) - score.getOpenPenalty(*(category+(lengthA-1)));
            } else {
                insert = _similarity_matrix[l-1][j] - score.getOpenPenalty(*(category+(lengthA-1)));
            }
            if (j==1 || _track_matrix[l][j - 1] == DELETION || _track_matrix[l][j - 1] == SNPORDELETION ||
                _track_matrix[l][j - 1] == SNPORINSERTIONORDELETION || _track_matrix[l][j-1] == INSERTIONORDELETION ) { //insertion
                del = _similarity_matrix[l][j - 1] - score.getExtendPenalty(*(category+(lengthA-1))) - score.getOpenPenalty(*(category+(lengthA-1)));
            } else {
                del = _similarity_matrix[l][j - 1] - score.getOpenPenalty(*(category+(lengthA-1)));
            }

            int32_t selected=0;
            if( del >insert && del==match  ){
                selected = del;
                _track_matrix[l][j] = SNPORDELETION;
            }else if( insert >del && insert == match  ){
                selected = match;
                _track_matrix[l][j] = SNPORINSERTION;
            }else if ( insert > match && insert > del){// prefer deletion
                int t = 1;
                while( l-t >=1 && (_track_matrix[l-t][j] == SNPORINSERTION || _track_matrix[l-t][j] == SNPORINSERTIONORDELETION || _track_matrix[l-t][j]==INSERTIONORDELETION ) ){
                    _track_matrix[l-t][j] = INSERTION;
                    ++t;
                }
                selected = insert;
                _track_matrix[l][j] = INSERTION;
            }else if( del > match && del > insert ){//prefer insertion, so that the INDELs could be put together
                int t = 1;
                while( j-t >=1 && (_track_matrix[l][j-t] == SNPORDELETION || _track_matrix[l][j-t] == SNPORINSERTIONORDELETION || _track_matrix[l][j-t]==INSERTIONORDELETION) ){
                    _track_matrix[l][j-t] = DELETION;
                    ++t;
                }
                selected = del;
                _track_matrix[l][j] = DELETION;
            }else if (match > insert && match > del){
                int t = 1;
                while( l-t >=1 && j-t>=1 && (_track_matrix[l-t][j-t] == SNPORINSERTION || _track_matrix[l-t][j-t] == SNPORINSERTIONORDELETION || _track_matrix[l-t][j-t]==SNPORDELETION ) ){
                    _track_matrix[l-t][j-t] = SNP;
                    ++t;
                }
                selected = match;
                _track_matrix[l][j] = SNP;
            }else if ( del >match && insert==del  ){
                selected = del;
                _track_matrix[l][j] = INSERTIONORDELETION;
            } else{
                selected = del;
                _track_matrix[l][j] = SNPORINSERTIONORDELETION;
            }
            _similarity_matrix[l][j] = selected;

        }
    }

    int32_t i = lengthB;
    int32_t j = lengthA;

    while (i > 0 || j > 0) {
        if (i == 0) {
            A.push('D');
            --j;
        } else if (j == 0) {
            A.push('I');
            --i;
        }else{
            if ( _track_matrix[i][j]==SNP || _track_matrix[i][j]==SNPORDELETION || _track_matrix[i][j]==SNPORINSERTION || _track_matrix[i][j]==SNPORINSERTIONORDELETION  ) {
                --i;
                --j;
                A.push('M');
            } else if (_track_matrix[i][j]==DELETION || _track_matrix[i][j]==INSERTIONORDELETION) {// Going to S(i, j-1) //deletion
                A.push('D');
                --j;
            } else {        //insertion
                A.push('I');
                --i;
            }
        }
    }


    for (i = 0; i <= lengthB; i++) {
        delete[] _similarity_matrix[i];
        delete[] _track_matrix[i];
    }
    delete[] _similarity_matrix;
    delete[] _track_matrix;
}
