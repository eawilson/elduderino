
#ifndef _MASH_H
#define _MASH_H

#include <sys/types.h>
#include <stdio.h>
#include <stdbool.h>

#ifndef _HASH_FUNCTION_T
#define _HASH_FUNCTION_T
typedef uint32_t (*hash_function_t)(const void *key, size_t length);
#endif


typedef struct mashentry_t {
    size_t key_size;
    size_t data_size;
    uint32_t next;
    } MashEntry;


typedef struct mashtable_t {
    uint32_t *buckets;
    uint32_t len_buckets;
    uint32_t buckets_occupied;
    MashEntry *entries;
    uint32_t len_entries;
    uint32_t entries_occupied;
    uint32_t *available_entries;
    void *keys;
    size_t max_key_len;
    void *data;
    size_t max_data_len;
    float bucket_resize;
    hash_function_t hash_function;
    } MashTable;



MashTable *mash_new(uint32_t size);
void mash_destroy(MashTable *ht);
int mash_put(MashTable *ht, const void *key, size_t key_size, const void *data, size_t data_size);
void *mash_get(MashTable *ht, const void *key, size_t key_size, size_t *data_size);
void *mash_pop(MashTable *ht, const void *key, size_t key_size, size_t *data_size);
void mash_summary_fprintf(MashTable *ht, FILE *fp);
void mash_contents_fprintf(MashTable *ht, FILE *fp);
void *mash_popall(MashTable *mt, const void **key, size_t *key_size, size_t *data_size, uint32_t *bucket);


#endif
