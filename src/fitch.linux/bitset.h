#include <stdint.h>

typedef int bitset_t;

/* Returns a pointer to a newly created bitset */
bitset_t * bitset_create(int max_size);

/* Frees memory allocated for bitset */
void bitset_destroy(bitset_t * s);

/* Adds element at given index */
int bitset_insert(bitset_t *s, int ind);

/* Removes element from given index */
int bitset_remove(bitset_t *s, int ind);

/* Sets all elements to 0 */
void bitset_clear(bitset_t *s);
 
/* Tests if given index is on or not */
int bitset_test(bitset_t *s, int ind);

/* Returns 1 if s is empty, else 0 */
int bitset_is_empty(bitset_t *s);

/* Returns number of elements in bitset */
int bitset_size(bitset_t *s);

/* Returns maximum possible elements in bitset */
int _bitset_max_size(bitset_t *s);

/* Returns number of words allocated for bitset */
int _bitset_words(bitset_t *s);

/* Performs bitwise OR to calculate set union */
void bitset_union(bitset_t *res, bitset_t *a, bitset_t *b);

/* Performs bitwise AND to calculate set intersection */
void bitset_intersection(bitset_t *res, bitset_t *a, bitset_t *b);

/* Prints all elements in the set */
void bitset_print(bitset_t *s);

int max(int a, int b);
int min(int a, int b);
