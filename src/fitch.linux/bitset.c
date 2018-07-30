#include "bitset.h"
#include <stdlib.h>
#include <stdio.h>
#define W_SIZE (8 * sizeof(int))

/* Returns pointer to a newly created bitset */
bitset_t * bitset_create(int max_size) {
    max_size --;
    bitset_t * temp = (bitset_t *)calloc((max_size / W_SIZE) + 2, sizeof(bitset_t));
    temp[0] = max_size / W_SIZE + 1;
    return temp;
}

/* Frees memory allocated for bitset */
void bitset_destroy(bitset_t * s) {
    free(s);
}

/* Adds element at given index */
int bitset_insert(bitset_t * s, int ind) {
    if (ind > _bitset_max_size(s)) return 1;
    s[(ind / W_SIZE) + 1] |= (1 << (ind % W_SIZE));
    return 0;
}

/* Removes element from given index */
int bitset_remove(bitset_t * s, int ind) {
    if (ind > _bitset_max_size(s)) return 1;
    s[(ind / W_SIZE) + 1] &= ~(1 << (ind % W_SIZE));
    return 0;
}

/* Sets all elements to 0 */
void bitset_clear(bitset_t * s) {
    for (int i=1; i < _bitset_words(s); i++) {
        s[i] &= 0;
    }
}
 
/* Tests if given index is on or not */
int bitset_test(bitset_t * s, int ind) {
    if (ind > _bitset_max_size(s)) return 0;
    return (s[(ind / W_SIZE) + 1] >> (ind % W_SIZE)) & 1;
}

/* Returns 1 if s is empty, else 0 */
int bitset_is_empty(bitset_t *s) {
    int n_words = s[0];
    int result = 0;
    for (int i = 1; i <= n_words; i++) {
        result |= s[i];
    }
    return result == 0;
}

/* Returns number of elements in bitset */
int bitset_size(bitset_t *s) {
    int size = 0;
    for (int i = 0; i < _bitset_max_size(s); i++) {
            size += bitset_test(s, i);
        }
    return size;
}

/* Returns maximum possible elements in bitset */
int _bitset_max_size(bitset_t *s) {
    return s[0] * W_SIZE;
}

/* Returns number of words allocated for bitset */
int _bitset_words(bitset_t *s) {
    return s[0];
}

/* Performs bitwise OR to calculate set union */
void bitset_union(bitset_t * res, bitset_t * a, bitset_t * b) {
    int a_words = _bitset_words(a);
    int b_words = _bitset_words(b);
    int min_words = min(a_words, b_words);
    for (int i = 1; i <= min_words; i++) {
        res[i] = a[i] | b[i];
    }
   if (a_words == b_words) return;
   if (a_words < b_words) {
       for (int i = min_words + 1; i < b_words; i++) res[i] = b[i];
    } else {
       for (int i = min_words + 1; i < a_words; i++) res[i] = a[i];
    }
}

/* Performs bitwise AND to calculate set intersection */
void bitset_intersection(bitset_t * res, bitset_t * a, bitset_t * b) {
    for (int i = 1; i <= min(_bitset_words(a), _bitset_words(b)); i++) {
        res[i] = a[i] & b[i];
    }
}

/* Prints all elements in the set */
void bitset_print(bitset_t *set) {
    for (int i = 0; i < _bitset_max_size(set); i++) {
        if (bitset_test(set, i)) printf("%d\n", i);
    }
}

int max(int a, int b) {
    return a > b ? a : b;
}

int min(int a, int b) {
    return a < b ? a : b;
}
