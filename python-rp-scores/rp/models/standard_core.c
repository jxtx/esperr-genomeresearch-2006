#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "standard_core.h"

#define bool int
#define false 0
#define true 1

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

/**
 * Return the flat size of a matrix having the given tuple size and radix
 */
int matrix_size( int order, int radix )
{
    return (int) pow( radix, order+1 );
}

/**
 * Return the index corresponding to the word string[offset-order:offset]
 */
int matrix_index( int order, int radix, int* string, int offset )
{
    int i, letter;

    int index = 0;
    int factor = 1;

    for ( i = 0; i < order + 1; i++ )
    {
        letter = string[ offset - i ];

        if ( letter < 0 || letter >= radix )
        {
            return -1;
        }

        index += letter * factor;
        factor *= radix;
    }

    return index;
}

int** new_counts( int order, int radix )
{
    int i, size;
    int** result = (int**) malloc( ( order + 1 ) * sizeof( int* ) );

    if ( result == NULL ) return NULL;

    for ( i = 0; i < order + 1; i++ )
    {
        size = matrix_size( i, radix ) * sizeof( int );
            
        result[i] = (int*) malloc( size );
        
        if ( result[i] == NULL ) return NULL;
        
        memset( result[i], 0, size );
    }

    return result;
}

void free_counts( int** v, int order )
{
    int i;

    for ( i = 0; i < order + 1; i++ )
    {
        free( v[i] );
    }

    free( v );
}

real* new_real_array( int size )
{
    real* result = ( real* ) malloc( size * sizeof( real ) );
    if ( result == NULL ) return NULL;
    memset( result, 0, size * sizeof( real ) );
    return result;
}

real** new_probs( int order, int radix )
{
    int i;
    real** result = (real**) malloc( ( order + 1 ) * sizeof( real* ) );
    if ( result == NULL ) return NULL;

    for ( i = 0; i < order + 1; i++ )
    {
        result[i] = new_real_array( matrix_size( i, radix ) );
    }

    return result;
}

void free_probs( real** v, int order)
{
    int i;

    for ( i = 0; i < order + 1; i++ )
    {
        free( v[i] );
    }

    free( v );
}


void fill_in_counts( int order, int radix, int** counts, int* string, int string_len )
{
    int i, j, index, max_lookback;

    int max_index = matrix_size( order, radix );

    for ( i = 0; i < string_len; i++ )
    {
        max_lookback = min( order, i ) + 1;
            
        for ( j = 0; j < max_lookback; j++ )
        {
            index = matrix_index( j, radix, string, i );
            
            if ( index == -1 )
            {
                continue;
            }
            
            assert( index >= 0 && index < max_index );
                    
            counts[ j ][ index ] += 1;
        }
    }
}

real** counts_to_probs( int order, int radix, int** counts, bool average, bool backoff )
{
    int i, row_start, size, prev_size, j, total, index;
    real prob;
    bool some_zero;

    real** probs = new_probs( order, radix );
    if ( probs == NULL ) return NULL;

    for ( i = 0; i < order+1; i++ )
    {
        size = matrix_size( i, radix );
        prev_size = matrix_size( i - 1, radix );

        for ( row_start = 0; row_start < size; row_start += radix )
        {
            total = 0;
            some_zero = false;

            for ( j = 0; j < radix; j++ )
            {
                index = row_start + j;

                assert( index >= 0 && index < size );

                if ( counts[ i ][ index ] == 0 )
                {
                    some_zero = true;
                }
                else
                {
                    total += counts[ i ][ index ];
                }
            }

            for ( j = 0; j < radix; j++ )
            {
                index = row_start + j;

                assert( index >= 0 && index < size );

                if ( backoff && total == 0 && i > 0 )
                {
                    prob = probs[ i - 1 ][ index % prev_size ];
                }
                else
                {
                    if ( some_zero )
                    {
                        prob = ( (real) counts[ i ][ index ] + 1 ) / ( (real) ( total + radix ) );
                    }
                    else
                    {
                        prob = ( (real) counts[ i ][ index ] ) / ( (real) total );
                    }
                }

		        if ( average && i > 0 )
		        {
		            probs[ i ][ index ] = prob + probs[ i - 1 ][ index % prev_size ];
		        }
		        else
		        {
		            probs[ i ][ index ] = prob;
		        }
            }
        }
    }

    return probs;
}

real* probs_to_score_matrix( int order, int radix, 
                             real** pos_probs, real** neg_probs, 
                             bool averaging )
{
    int i;

    int size = matrix_size( order, radix );
    real* scores = new_real_array( size );
    real log_order = log( order );

    if ( scores == NULL ) return NULL;

    for ( i = 0; i < size; i++ )
    {
        // If averaging the probs need to be divided by the order
        
        if ( averaging )
        {
            scores[i] = ( log( pos_probs[order][i] ) - log_order ) 
                      - ( log( neg_probs[order][i] ) - log_order );   
        }
        else
        {
            scores[i] = log( pos_probs[order][i] ) 
                      - log( neg_probs[order][i] );
        }
    }

    return scores;
}

void free_scores( real* scores )
{
    free( scores );
}

bool score_string( int order, int radix, real* score_matrix, 
                   int* string, int start, int length, real* rval )
{
    int i, index;
    real score = 0;
    int valid_tuples = 0;
    int stop = start + length;

    for ( i = start + order; i < stop; i++ )
    {
        index = matrix_index( order, radix, string, i );

        // Skip tuples containing invalid symbols
        if ( index == -1 )
        {
            continue;
        }
        else 
        {
            valid_tuples++;
        }

        score += score_matrix[ index ];
    }

    if ( valid_tuples > 0 )
    {
        *rval = score / (real) valid_tuples;
	    return true;
    }
    else
    {
        *rval = 0;
        return false;
    }
}

bool score_string_positions( int order, int radix, real* score_matrix, 
                             int* string, float* target, int start, int length )
{
    int i, index;
    real score = 0;
    int valid_tuples = 0;
    int stop = start + length;

    for ( i = start + order; i < stop; i++ )
    {
        index = matrix_index( order, radix, string, i );

        // Skip tuples containing invalid symbols
        if ( index == -1 )
        {
            continue;
        }
        else 
        {
            valid_tuples++;
        }

        target[i] = score_matrix[ index ];
    }
}
