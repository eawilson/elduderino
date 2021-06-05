
#ifndef _HASH_H
#define _HASH_H

#include <sys/types.h>
#include <stdio.h>
#include <stdbool.h>


#ifndef _HASH_FUNCTION_T
#define _HASH_FUNCTION_T
typedef uint32_t (*hash_function_t)(const void *key, size_t length);
#endif


typedef struct hashentry_t {
    const void *key;
    size_t key_size;
    const void *data;
    size_t data_size;
    uint32_t next;
    } HashEntry;


typedef struct hashtable_t {
    uint32_t *buckets;
    uint32_t len_buckets;
    uint32_t buckets_occupied;
    HashEntry *entries;
    uint32_t len_entries;
    uint32_t entries_occupied;
    uint32_t *available_entries;
    float bucket_resize;
    hash_function_t hash_function;
    } HashTable;



HashTable *hash_new(uint32_t size);
void hash_destroy(HashTable *ht);
int hash_put(HashTable *ht, const void *key, size_t key_size, const void *data, size_t data_size);
void *hash_get(HashTable *ht, const void *key, size_t key_size, size_t *data_size);
void *hash_pop(HashTable *ht, const void *key, size_t key_size, size_t *data_size);
void hash_summary_fprintf(HashTable *ht, FILE *fp);
void hash_contents_fprintf(HashTable *ht, FILE *fp);
int hash_validate(HashTable *ht, FILE *fp);


#endif
