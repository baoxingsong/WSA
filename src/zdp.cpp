
#include "zdp.h"


/*
 * =====================================================================================
 *
 *       Filename:  zdp.cpp
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


a weighted dynamic programming sequence alignment method ZDP (Zebric dynamic programming)

The standard needleman wunsch algorithm set up a score matric and generate alignment by tracing back
 Since we have different score stragies for each basepair the tracing back step would query the weighted score again, which is time cosuming
 So here when set up the score matrix, we setup a trace matrix at the same time, and we coould trace back using the trace matrix

 ************************************************************************/




void zdp(const char * seqAChar, const char * seqBChar, const int32_t & lengthA, const int32_t & lengthB, const int16_t * category, std::stack<char> & A, Score & score){
    if( lengthA==0 || 0==lengthB ){ // if the lengght any of the reference sequence or query sequence is 0, we do not need to align them.
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

    // internalize the score matrix begin
    auto ** _similarity_matrix = new int32_t*[lengthB + 1];
    for (int i = 0; i < (lengthB + 1); i++) {
        _similarity_matrix[i] = new int32_t[lengthA + 1];
    }
    for (int i = 0; i<= lengthB; i++) {
        for (int j = 0; j<= 1; j++) {
            _similarity_matrix[i][j] = - score.getExtendPenalty(*(category)) - score.getOpenPenalty(*(category));
        }
    }
    for (int j = 0; j<= lengthA; j++){
        for (int i = 0; i<= 1; i++) {
            _similarity_matrix[i][j] = - score.getExtendPenalty(*(category)) - score.getOpenPenalty(*(category));
        }
    }
    _similarity_matrix[0][0] = 0;
    // internalize the score matrix end

    // internalize the _trace_matrix matric begin
    auto **_trace_matrix = new VARIANTCATEGORY*[lengthB + 1];
    for (int i = 0; i < (lengthB + 1); i++) {
        _trace_matrix[i] = new VARIANTCATEGORY[lengthA + 1];
    }

    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= lengthA; j++) {
            _trace_matrix[i][j] = INSERTION;
        }
    }
    for (int j = 0; j <= 1; j++){
        for (int i = 0; i <= lengthB; i++) {
            _trace_matrix[i][j] = DELETION;
        }
    }
    _trace_matrix[0][0] = SNPORINSERTIONORDELETION; // we should never use this matrix cell
    // internalize the _trace_matrix matric end

    int l;
    int match = 0, insert = 0, del = 0;
    for (int j=1; j < lengthA + 1; ++j){
        for (l = 1; l < lengthB + 1; ++l) {
            match = _similarity_matrix[l - 1][j - 1] + score.getScore(*(category+(j-1)), *(seqAChar+j-1), *(seqBChar+l-1));

            if (_trace_matrix[l-1][j] == INSERTION || _trace_matrix[l-1][j] == SNPORINSERTION
                || _trace_matrix[l-1][j] == SNPORINSERTIONORDELETION || _trace_matrix[l - 1][j] == INSERTIONORDELETION ) { //deletion
                insert = _similarity_matrix[l-1][j] - score.getExtendPenalty(*(category+(j-1))) - score.getOpenPenalty(*(category+(j-1)));
            } else {
                insert = _similarity_matrix[l-1][j] - score.getOpenPenalty(*(category+(j-1)));
            }
            if ( _trace_matrix[l][j - 1] == DELETION || _trace_matrix[l][j - 1] == SNPORDELETION ||
                _trace_matrix[l][j - 1] == SNPORINSERTIONORDELETION || _trace_matrix[l][j-1] == INSERTIONORDELETION ) { //insertion
                del = _similarity_matrix[l][j - 1] - score.getExtendPenalty(*(category+(j-1))) - score.getOpenPenalty(*(category+(j-1)));
            } else {
                del = _similarity_matrix[l][j - 1] - score.getOpenPenalty(*(category+(j-1)));
            }

            // check which of SNP/INSERTION/DELETION is largest and fill the score matric and trace matric
            int32_t selected=0;
            if( del >insert && del==match  ){
                selected = del;
                _trace_matrix[l][j] = SNPORDELETION;
            }else if( insert >del && insert == match  ){
                selected = match;
                _trace_matrix[l][j] = SNPORINSERTION;
            }else if ( insert > match && insert > del){// prefer deletion
                int t = 1;
                // if for this cell the insertion is the largest, the previouse cell should be updated as insertion
                while( l-t >=1 && (_trace_matrix[l-t][j] == SNPORINSERTION || _trace_matrix[l-t][j] == SNPORINSERTIONORDELETION || _trace_matrix[l-t][j]==INSERTIONORDELETION ) ){
                    _trace_matrix[l-t][j] = INSERTION;
                    ++t;
                }
                selected = insert;
                _trace_matrix[l][j] = INSERTION;
            }else if( del > match && del > insert ){//prefer insertion, so that the INDELs could be put together
                int t = 1;
                while( j-t >=1 && (_trace_matrix[l][j-t] == SNPORDELETION || _trace_matrix[l][j-t] == SNPORINSERTIONORDELETION || _trace_matrix[l][j-t]==INSERTIONORDELETION) ){
                    _trace_matrix[l][j-t] = DELETION;
                    ++t;
                }
                selected = del;
                _trace_matrix[l][j] = DELETION;
            }else if (match > insert && match > del){
                int t = 1;
                while( l-t >=1 && j-t>=1 && (_trace_matrix[l-t][j-t] == SNPORINSERTION || _trace_matrix[l-t][j-t] == SNPORINSERTIONORDELETION || _trace_matrix[l-t][j-t]==SNPORDELETION ) ){
                    _trace_matrix[l-t][j-t] = SNP;
                    ++t;
                }
                selected = match;
                _trace_matrix[l][j] = SNP;
            }else if ( del >match && insert==del  ){
                selected = del;
                _trace_matrix[l][j] = INSERTIONORDELETION;
            } else{
                selected = del;
                _trace_matrix[l][j] = SNPORINSERTIONORDELETION;
            }
            _similarity_matrix[l][j] = selected;

        }
    }

    // begin to trace back using the _trace_matrix
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
            if ( _trace_matrix[i][j]==SNP || _trace_matrix[i][j]==SNPORDELETION || _trace_matrix[i][j]==SNPORINSERTION || _trace_matrix[i][j]==SNPORINSERTIONORDELETION  ) {
                --i;
                --j;
                A.push('M');
            } else if (_trace_matrix[i][j]==DELETION || _trace_matrix[i][j]==INSERTIONORDELETION) {// Going to S(i, j-1) //deletion
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
        delete[] _trace_matrix[i];
    }
    delete[] _similarity_matrix;
    delete[] _trace_matrix;
}
