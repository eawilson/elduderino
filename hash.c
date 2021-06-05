#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

#include "hash.h"


// static uint32_t perl_hash(const void *key, size_t length);
static uint32_t fnv1a_hash(const void *key, size_t length);
static int resize_buckets(HashTable *ht, uint32_t size);
static int resize_entries(HashTable *ht, uint32_t size);



static uint32_t fnv1a_hash(const void *key, size_t length) {
    register size_t i = 0;
    register uint32_t hash = 2166136261;
    register const uint32_t fnv_prime = 16777619;
    for (i = 0; i < length; ++i) {
        hash = (hash ^ ((unsigned char *)key)[i]) * fnv_prime;
        }
    
    return hash;
    }



static int cmp_uint32(const void *p1, const void *p2) {
    return *(uint32_t *)p1 - *(uint32_t *)p2;
    }



// static uint32_t perl_hash(const void *key, size_t length) {
//     register size_t i = length;
//     register uint32_t hash = 0;
//     register const unsigned char *s = (unsigned char *)key;
//     
//     while (i--) {
//         hash += *s++;
//         hash += (hash << 10);
//         hash ^= (hash >> 6);
//         }
//     hash += (hash << 3);
//     hash ^= (hash >> 11);
//     hash += (hash << 15);
// 
//     return hash;
//     }



HashTable *hash_new(uint32_t size) {
    HashTable *ht = NULL;
    
    size = size > 0 ? size : 1;
    
    if ((ht = (HashTable *)calloc(1, sizeof(HashTable))) == NULL) {
        return NULL;
        }
    
    ht->bucket_resize = 0.7;
    ht->hash_function = &fnv1a_hash;
    
    if (resize_buckets(ht, size) == -1 || resize_entries(ht, ht->len_buckets) == -1) {
        hash_destroy(ht);
        return NULL;
        }
    return ht;
    }



void hash_destroy(HashTable *ht) {
    free(ht->buckets);
    free(ht->entries);
    free(ht->available_entries);
    free(ht);
    }



void hash_summary_fprintf(HashTable *ht, FILE *fp) {
    fprintf(fp, "%li / %li Buckets\n", (long)ht->buckets_occupied, (long)ht->len_buckets);
    fprintf(fp, "%li / %li Entries\n", (long)ht->entries_occupied, (long)ht->len_entries);
    }



void hash_contents_fprintf(HashTable *ht, FILE *fp) {
    uint32_t bucket = 0, i = 0;
    HashEntry *entry = NULL;
    
    for (bucket = 0; bucket < ht->len_buckets; ++bucket) {
        if ((i = ht->buckets[bucket]) != UINT32_MAX) {
            fprintf(fp, "Bucket: %li\n", (long)bucket);
            for (; i != UINT32_MAX; i = entry->next) {
                entry = ht->entries + i;
                fprintf(fp, "    %.*s    %.*s\n", (int)entry->key_size, (char *)entry->key, (int)entry->data_size, (char *)entry->data);
                }
            }
        }
    }



int hash_put(HashTable *ht, const void *key, size_t key_size, const void *data, size_t data_size) {
    uint32_t bucket = 0, i = 0;
    HashEntry *entry = NULL;
    
    //fprintf(stderr, "put\n");
    if (ht->buckets_occupied > ht->len_buckets * ht->bucket_resize) {
        if ((resize_buckets(ht, ht->len_buckets * 2)) == -1) {
            return -1;
            }
        }
    
    if (ht->entries_occupied == ht->len_entries) {
        if ((resize_entries(ht, ht->len_entries * 2)) == -1) {
            return -1;
            }
        }
    
    bucket = ht->hash_function(key, key_size) & (ht->len_buckets - 1);
    if (ht->buckets[bucket] == UINT32_MAX) {
        ++ht->buckets_occupied;
        }
    
    i = ht->available_entries[ht->entries_occupied++];
    entry = ht->entries + i;
    entry->key = key;
    entry->key_size = key_size;
    entry->data = data;
    entry->data_size = data_size;
    entry->next = ht->buckets[bucket];
    ht->buckets[bucket] = i;
    return 0;
    }



void *hash_get(HashTable *ht, const void *key, size_t key_size, size_t *data_size) {
    uint32_t i = 0, bucket = 0;
    HashEntry *entry = NULL;
    
    //fprintf(stderr, "get\n");
    bucket = ht->hash_function(key, key_size) & (ht->len_buckets - 1);
    for (i = ht->buckets[bucket]; i != UINT32_MAX; i = entry->next) {
        entry = ht->entries + i;
        if (key_size == entry->key_size && memcmp(key, entry->key, key_size) == 0) {
            *data_size = entry->data_size;
            return (void *)entry->data;
            }
        }
    return NULL;
    }



void *hash_pop(HashTable *ht, const void *key, size_t key_size, size_t *data_size) {
    uint32_t i = 0, bucket = 0;
    HashEntry *entry = NULL, *previous = NULL;
    
    //fprintf(stderr, "pop\n");
    bucket = ht->hash_function(key, key_size) & (ht->len_buckets - 1);
    for (i = ht->buckets[bucket]; i != UINT32_MAX; i = entry->next) {
        entry = ht->entries + i;
        if (key_size == entry->key_size && memcmp(key, entry->key, key_size) == 0) {
            if (previous == NULL) {
                if ((ht->buckets[bucket] = entry->next) == UINT32_MAX) {
                    --ht->buckets_occupied;
                    }
                }
            else {
                previous->next = entry->next;
                }
            
            ht->available_entries[--ht->entries_occupied] = i;
            
            *data_size = entry->data_size;
            void *data = (void *)entry->data;
            memset(entry, 0, sizeof(HashEntry)); // this could be changed to just zero *key in production code but_vaidate would fail in this situation
            return data;
            }
        
        previous = entry;
        }
    return NULL;
    }



int hash_validate(HashTable *ht, FILE *fp) {
    uint32_t i = 0, buckets_occupied = 0, *available = NULL, j = 0;
    HashEntry *entry = NULL, zero_entry;
    int ret = 0;
    memset(&zero_entry, 0, sizeof(HashEntry));
   
    // check buckets_occupied <= len_buckets
    if (ht->buckets_occupied > ht->len_buckets) {
        fprintf(fp, "buckets occupied (%u) > len buckets (%u)\n", (unsigned)ht->buckets_occupied, (unsigned)ht->len_buckets);
        ret = -1;
        }
    
    
    // check buckets_occupied is correct and that only valid entries are referenced from the buckets
    for (i = 0; i < ht->len_buckets; ++i) {
        if (ht->buckets[i] < UINT32_MAX) {
            ++buckets_occupied;
            if (ht->buckets[i] >= ht->len_entries) {
                fprintf(fp, "entry referenced in bucket (%u) >= len entries (%u)\n", (unsigned)ht->buckets[i], (unsigned)ht->len_entries);
                ret = -1;
                break;
                }
            else if (memcmp(&zero_entry, ht->entries + ht->buckets[i], sizeof(HashEntry)) == 0) {
                fprintf(fp, "zeroed entry referenced in bucket\n");
                ret = -1;
                break;
                }
            }
        }
    if (buckets_occupied != ht->buckets_occupied) {
        fprintf(fp, "actual buckets occupied (%u) != buckets occupied (%u)\n", (unsigned)buckets_occupied, (unsigned)ht->buckets_occupied);
        ret = -1;
        }
    
    
    // check entries_occupied <= len_entries
    if (ht->entries_occupied > ht->len_entries) {
        fprintf(fp, "entries occupied (%u) > len entries (%u)\n", (unsigned)ht->entries_occupied, (unsigned)ht->len_entries);
        ret = -1;
        }
    
    
    // create an array of sorted available entries for next checks
    if ((available = (uint32_t *)malloc((ht->len_entries - ht->entries_occupied) * sizeof(uint32_t))) == NULL) {
        fprintf(fp, "FAILED TO ALLOCATE MEMORY TO COMPLETE VALIDATION\n");
        return -2;
        }
    memcpy(available, ht->available_entries + ht->entries_occupied, (ht->len_entries - ht->entries_occupied) * sizeof(uint32_t));
    qsort(available, ht->len_entries - ht->entries_occupied, sizeof(uint32_t), cmp_uint32);
    
    
    // check no duplicate available entries
    for (i = 0; i < ht->len_entries - ht->entries_occupied - 1; ++i) {
        if (available[i] == available[i + 1]) {
            fprintf(fp, "duplicate available entry %u\n", (unsigned)available[i]);
            ret = -1;
            break;
            }
        }
    
    
    j = 0; // index of next available
    for (i = 0; i < ht->len_entries; ++i) {
        // available entry therefore check if zeroed
        if (j < ht->len_entries - ht->entries_occupied && i == available[j]) {
            if (memcmp(&zero_entry, ht->entries + i, sizeof(HashEntry) != 0)) {
                fprintf(fp, "unused entry not zeroed\n");
                ret = -1;
                break;
                }
            ++j;
            }
        // used entry therefore check valid
        else {
            entry = ht->entries +i;
            if (memcmp(&zero_entry, entry, sizeof(HashEntry) == 0)) {
                fprintf(fp, "used entry zeroed\n");
                ret = -1;
                break;
                }
            else if (entry->key == NULL || entry->data == NULL || entry->key_size == 0 || entry->data_size == 0) {
                fprintf(fp, "entry has null pointers or zero length elements\n");
                ret = -1;
                break;
                }
            else if (entry->next != UINT32_MAX && memcmp(&zero_entry, ht->entries + entry->next, sizeof(HashEntry) == 0)) {
                fprintf(fp, "entry next references a zeroed element\n");
                ret = -1;
                break;
                }
            }
        }
    if (j < ht->len_entries - ht->entries_occupied) {
        fprintf(fp, "invalid available entries\n");
        ret = -1;
        }
        
        
    if (ret == 0) {
        fprintf(fp, "valid\n");
        }
    
    
    free(available);
    return ret;
    }


static int resize_buckets(HashTable *ht, uint32_t size) {
    uint32_t bucket = 0, i = 0;
    HashEntry *entry = NULL;
    
    //fprintf(stderr, "resize buckets\n");
    ht->len_buckets = 1;
    while (ht->len_buckets < size) {
        ht->len_buckets <<= 1;
        }
    free(ht->buckets);
    if ((ht->buckets = (uint32_t *)malloc(ht->len_buckets * sizeof(uint32_t))) == NULL) {
        return -1;
        }
    memset(ht->buckets, 0xFF, ht->len_buckets * sizeof(uint32_t));
    
    ht->buckets_occupied = 0;
    for (i = 0; i < ht->len_entries; ++i) {
        entry = ht->entries + i;
        if (entry->key != NULL) {
            bucket = ht->hash_function(entry->key, entry->key_size) & (ht->len_buckets - 1);
            if (ht->buckets[bucket] == UINT32_MAX) {
                ++ht->buckets_occupied;
                }
            entry->next = ht->buckets[bucket];
            ht->buckets[bucket] = i;
            }
        }
    return 0;
    }



static int resize_entries(HashTable *ht, uint32_t size) {
    uint32_t i = 0;
    void *ptr = NULL;
    
    if (size > ht->len_entries) {
        // fprintf(stderr, "resize entries\n");
        if ((ptr = realloc(ht->entries, size * sizeof(HashEntry))) == NULL) {
            return -1;
            }
        ht->entries = (HashEntry *)ptr;
        memset(ht->entries + ht->len_entries, 0, (size - ht->len_entries) * sizeof(HashEntry));
        
        if ((ptr = realloc(ht->available_entries, size * sizeof(uint32_t))) == NULL) {
            return -1;
            }
        ht->available_entries = (uint32_t *)ptr;
        for (i = ht->len_entries; i < size; ++i) {
            ht->available_entries[i] = i;
            }
        
        ht->len_entries = size;
        }
        
    return 0;
    }


