
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "simple_periodic_core.h"

#define bool int
#define false 0
#define true 1

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

#undef NDEBUG
#include <assert.h>

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

int* new_counts( int order, int radix )
{
    int size = matrix_size( order, radix ) * sizeof( int );
    int * result = (int*) malloc( size );
    if ( result == NULL ) return NULL;
    memset( result, 0, size );
    return result;
}

int** new_counts_2( int order1, int order2, int radix1, int radix2 )
{
    int i;
    int nmodels = matrix_size( order2, radix1 );
    int size = matrix_size( order2, radix2 ) * sizeof( int );
    
    int ** result = (int**) malloc( nmodels * sizeof( int * ) );
    if ( result == NULL ) return NULL;
    
    for ( i = 0; i < nmodels; i++ )
    {
        result[i] = (int*) malloc( size );
        if ( result[i] == NULL ) return NULL;
        memset( result[i], 0, size );
    }
    
    // fprintf( stderr, "%d, %d\n", nmodels, size/sizeof(int) );
    // fflush( stderr ); 
    
    return result;
}

void free_counts( int* v, int order )
{
    free( v );
}

real* new_real_array( int size )
{
    real* result = ( real* ) malloc( size * sizeof( real ) );
    if ( result == NULL ) return NULL;
    memset( result, 0, size * sizeof( real ) );
    return result;
}

real* new_probs( int order, int radix )
{

    return new_real_array( matrix_size( order, radix ) );
}

void free_probs( real* v, int orde)
{
    free( v );
}

/**
 * periodic models: only positions where "( index + phase ) % period" will 
 * be counted
 *
 * phase = 0, period = 1 gives us the normal model
 */
void fill_in_counts( int order, int radix, int* counts, int* string, int string_len, int phase, int period )
{
    int i, index;

    int max_index = matrix_size( order, radix );

    for ( i = phase; i < string_len; i += period )
    {
        if ( i < order )
        {
            continue;
        }
        
        index = matrix_index( order, radix, string, i );
            
        if ( index == -1 )
        {
            continue;
        }
            
        assert( index >= 0 && index < max_index );
                    
        counts[ index ] += 1;
    }
}

void fill_in_counts_2( int order1, int order2, int radix1, int radix2, 
                       int** counts, int* string1, int* string2,
                       int string_len, int phase, int period )
{
    int i, model_index, transition_index;

    int max_index1 = matrix_size( order1, radix1 );
    int max_index2 = matrix_size( order2, radix2 );

    for ( i = phase; i < string_len; i += period )
    {
        if ( i < order1 || i < order2 )
        {
            continue;
        }
        
        model_index = matrix_index( order2, radix1, string1, i );
        
        transition_index = matrix_index( order2, radix2, string2, i );   
            
        if ( model_index == -1 || transition_index == -1 )
        {
            continue;
        }
            
        assert( model_index >= 0 && model_index < max_index1 );
        assert( transition_index >= 0 && transition_index < max_index2 );
        
        // fprintf( stderr, "fill: %d, %d, %d, %d\n",  order1,  order2,  radix1,  radix2 );   
        // fprintf( stderr, "----: %d, %d\n", model_index, transition_index );
        // fflush( stderr );    
                    
        counts[ model_index ][ transition_index ] += 1;
    }
}

real* counts_to_probs( int order, int radix, int* counts )
{
    int j, row_start, size, total, index;
    real prob;
    bool some_zero;

    real* probs = new_probs( order, radix );
    if ( probs == NULL ) return NULL;

    size = matrix_size( order, radix );

    for ( row_start = 0; row_start < size; row_start += radix )
    {
        total = 0;
        some_zero = false;

        for ( j = 0; j < radix; j++ )
        {
            index = row_start + j;

            assert( index >= 0 && index < size );

            if ( counts[ index ] == 0 )
            {
                some_zero = true;
            }
            else
            {
                total += counts[ index ];
            }
        }

        for ( j = 0; j < radix; j++ )
        {
            index = row_start + j;

            assert( index >= 0 && index < size );

            if ( some_zero )
            {
                prob = ( (real) counts[ index ] + 1 ) / ( (real) ( total + radix ) );
            }
            else
            {
                prob = ( (real) counts[ index ] ) / ( (real) total );
            }

	        probs[ index ] = prob;
        }
    }

    return probs;
}

real* probs_to_score_matrix( int order, int radix, 
                             real* pos_probs, real* neg_probs )
{
    int i;

    int size = matrix_size( order, radix );
    real* scores = new_real_array( size );

    if ( scores == NULL ) return NULL;

    for ( i = 0; i < size; i++ )
    {
        scores[i] = log( pos_probs[i] ) - log( neg_probs[i] );
    }

    return scores;
}

void free_scores( real* scores )
{
    free( scores );
}

bool score_string( int order, int radix, real** score_matrix, 
                   int* string, int start, int length, real* rval, int period )
{
    int i, index, phase;
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
        
        phase = i % period;

        score += score_matrix[ phase ][ index ];
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

bool score_string_2( int order1, int order2, int radix1, int radix2, 
                     real** score_matrix1, real *** score_matrix2,
                     int* string1, int* string2, 
                     int start, int length, real* rval, int period )
{
    int i, index1, index1b, index2, phase;
    real score = 0;
    int valid_tuples = 0;
    int stop = start + length;

    for ( i = start + order1; i < stop; i++ )
    {
        index1 = matrix_index( order1, radix1, string1, i );
        index1b = matrix_index( order2, radix1, string1, i );
        index2 = matrix_index( order2, radix2, string2, i );

        // fprintf( stderr, "foo: %d, %d\n", index1, index2 ); fflush( stderr );

        // Skip tuples containing invalid symbols
        if ( index1 == -1 || index2 == -1 )
        {
            continue;
        }
        else 
        {
            valid_tuples++;
        }
        
        phase = i % period;

        score += score_matrix1[ phase ][ index1 ];
        score += score_matrix2[ phase ][ index1b ][ index2 ];
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

/*
bool score_string_positions( int order, int radix, real* score_matrix, 
                             int* string, float* target, int start, int length )
{
    int i, index;
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
    
    return true;
}
*/