#ifndef __TRAIN_HELPER_H__
#define __TRAIN_HELPER_H__

typedef double real;

int** new_counts( int tsize, int radix );
void free_counts( int** v, int tsize );
void fill_in_counts( int tsize, int radix, int** counts, int* string, int string_len );
real** counts_to_probs( int order, int radix, int** counts );
void free_probs( real** probs, int order );
real* new_real_array( int size );
real* probs_to_score_matrix( int tsize, int radix, real** pos_probs, real** neg_probs );
void free_scores( real* scores );

#endif
