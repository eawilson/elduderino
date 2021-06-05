#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>

#include "mash.h"



// static uint32_t perl_hash(const void *key, size_t length);
static uint32_t fnv1a_hash(const void *key, size_t length);
static int resize_buckets(MashTable *mt, uint32_t size);
static int resize_entries(MashTable *mt, uint32_t size);
static int resize_keys(MashTable *mt, size_t max_key_len);
static int resize_data(MashTable *mt, size_t max_data_len);



static uint32_t fnv1a_hash(const void *key, size_t length) {
    register size_t i = 0;
    register uint32_t hash = 2166136261;
    register const uint32_t fnv_prime = 16777619;
    
    for (i = 0; i < length; ++i) {
        hash = (hash ^ ((unsigned char *)key)[i]) * fnv_prime;
        }
    
    return hash;
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



MashTable *mash_new(uint32_t size) {
    MashTable *mt = NULL;
    
    if ((mt = (MashTable *)calloc(1, sizeof(MashTable))) == NULL) {
        return NULL;
        }
    
    mt->bucket_resize = 0.7;
    mt->hash_function = &fnv1a_hash;
    
    if (resize_buckets(mt, size) == -1 || resize_entries(mt, mt->len_buckets) == -1) {
        mash_destroy(mt);
        return NULL;
        }
    return mt;
    }



void mash_destroy(MashTable *mt) {
    free(mt->buckets);
    free(mt->entries);
    free(mt->available_entries);
    free(mt->keys);
    free(mt->data);
    free(mt);
    }



void mash_summary_fprintf(MashTable *mt, FILE *fp) {
    fprintf(fp, "%li / %li Buckets\n", (long)mt->buckets_occupied, (long)mt->len_buckets);
    fprintf(fp, "%li / %li Entries\n", (long)mt->entries_occupied, (long)mt->len_entries);
    }



void mash_contents_fprintf(MashTable *mt, FILE *fp) {
    uint32_t bucket = 0, i = 0;
    MashEntry *entry = NULL;
    
    for (bucket = 0; bucket < mt->len_buckets; ++bucket) {
        if ((i = mt->buckets[bucket]) != UINT32_MAX) {
            fprintf(fp, "Bucket: %li\n", (long)bucket);
            for (; i != UINT32_MAX; i = entry->next) {
                entry = mt->entries + i;
                fprintf(fp, "    %.*s    %.*s\n", (int)entry->key_size, (char *)(mt->keys + (i * mt->max_key_len)),
                                                  (int)entry->data_size, (char *)(mt->data + (i * mt->max_data_len)));
                }
            }
        }
    }



int mash_put(MashTable *mt, const void *key, size_t key_size, const void *data, size_t data_size) {
    uint32_t bucket = 0, i = 0;
    MashEntry *entry = NULL;
        
    // Key size 0 is used to identify empty buckets
    if (key_size == 0) {
        return -1;
        }
    
    if (mt->buckets_occupied > mt->len_buckets * mt->bucket_resize) {
        if ((resize_buckets(mt, mt->len_buckets * 2)) == -1) {
            return -1;
            }
        }
    
    if (mt->entries_occupied == mt->len_entries) {
        if ((resize_entries(mt, mt->len_entries * 2)) == -1) {
            return -1;
            }
        }
    
    if (key_size > mt->max_key_len) {
        if ((resize_keys(mt, key_size)) == -1) {
            return -1;
            }
        }
    
    if (data_size > mt->max_data_len) {
        if ((resize_data(mt, data_size)) == -1) {
            return -1;
            }
        }
    
    bucket = mt->hash_function(key, key_size) & (mt->len_buckets - 1);
    if (mt->buckets[bucket] == UINT32_MAX) {
        ++mt->buckets_occupied;
        }
    
    i = mt->available_entries[mt->entries_occupied++];
    memcpy(mt->keys + (i * mt->max_key_len), key, key_size);
    memcpy(mt->data + (i * mt->max_data_len), data, data_size);
    entry = mt->entries + i;
    entry->key_size = key_size;
    entry->data_size = data_size;
    entry->next = mt->buckets[bucket];
    mt->buckets[bucket] = i;
    return 0;
    }



void *mash_get(MashTable *mt, const void *key, size_t key_size, size_t *data_size) {
    uint32_t i = 0, bucket = 0;
    MashEntry *entry = NULL;
    
    bucket = mt->hash_function(key, key_size) & (mt->len_buckets - 1);
    i = mt->buckets[bucket];
    while (i != UINT32_MAX) {
        entry = mt->entries + i;
        if (key_size == entry->key_size && memcmp(key, mt->keys + (i * mt->max_key_len), key_size) == 0) {
            *data_size = entry->data_size;
            return mt->data + (i * mt->max_data_len);
            }
        i = entry->next;
        }
    return NULL;
    }



void *mash_pop(MashTable *mt, const void *key, size_t key_size, size_t *data_size) {
    uint32_t i = 0, bucket = 0;
    MashEntry *entry = NULL, *previous = NULL;
    
    bucket = mt->hash_function(key, key_size) & (mt->len_buckets - 1);
    i = mt->buckets[bucket];
    while (i != UINT32_MAX) {
        entry = mt->entries + i;
        if (key_size == entry->key_size && memcmp(key, mt->keys + (i * mt->max_key_len), key_size) == 0) {
            if (previous == NULL) {
                if ((mt->buckets[bucket] = entry->next) == UINT32_MAX) {
                    --mt->buckets_occupied;
                    }
                }
            else {
                previous->next = entry->next;
                }
            
            mt->available_entries[--mt->entries_occupied] = i;
            entry->key_size = 0;
            *data_size = entry->data_size;
            return mt->data + (i * mt->max_data_len);
            }
        
        previous = entry;
        i = entry->next;
        }
    return NULL;
    }



void *mash_popall(MashTable *mt, const void **key, size_t *key_size, size_t *data_size, uint32_t *bucket) {
    uint32_t i = 0;
    MashEntry *entry = NULL, *previous = NULL;
    
    for (; *bucket < mt->len_buckets; ++*bucket) {
        while ((i = mt->buckets[*bucket]) != UINT32_MAX) {
            for(entry = NULL; i != UINT32_MAX; i = entry->next) {
                previous = entry;
                entry = mt->entries + i;
                if (*key == NULL) {
                    *key = mt->keys + (i * mt->max_key_len);
                    *key_size = entry->key_size;
                    }
                else if (*key_size != entry->key_size || memcmp(*key, mt->keys + (i * mt->max_key_len), *key_size) != 0) {
                    continue;
                    }
                
                if (previous == NULL) {
                    if ((mt->buckets[*bucket] = entry->next) == UINT32_MAX) {
                        --mt->buckets_occupied;
                        }
                    }
                else {
                    previous->next = entry->next;
                    }
                
                mt->available_entries[--mt->entries_occupied] = i;
                entry->key_size = 0;
                
                *data_size = entry->data_size;
                return mt->data + (i * mt->max_data_len);

                }
            
            *key = NULL;
            }
        }
    return NULL;
    }



static int resize_buckets(MashTable *mt, uint32_t size) {
    uint32_t bucket = 0, i = 0;
    MashEntry *entry = NULL;
    
    mt->len_buckets = 1;
    while (mt->len_buckets < size) {
        mt->len_buckets <<= 1;
        }

    free(mt->buckets);
    if ((mt->buckets = (uint32_t *)malloc(mt->len_buckets * sizeof(uint32_t))) == NULL) {
        return -1;
        }
    memset(mt->buckets, 0xFF, mt->len_buckets * sizeof(uint32_t));

    if (mt->buckets_occupied) {
        mt->buckets_occupied = 0;
        for (i = 0; i < mt->len_entries; ++i) {
            entry = mt->entries + i;
            if (entry->key_size != 0) {
                bucket = mt->hash_function(mt->keys + (i * mt->max_key_len), entry->key_size) & (mt->len_buckets - 1);
                if (mt->buckets[bucket] == UINT32_MAX) {
                    ++mt->buckets_occupied;
                    }
                entry->next = mt->buckets[bucket];
                mt->buckets[bucket] = i;
                }
            }
        }
    
    return 0;
    }



static int resize_entries(MashTable *mt, uint32_t size) {
    uint32_t i = 0;
    void *ptr = NULL;
    
    if (size > mt->len_entries) {
        if ((ptr = realloc(mt->entries, size * sizeof(MashEntry))) == NULL) {
            return -1;
            }
        mt->entries = (MashEntry *)ptr;
        memset(mt->entries + mt->len_entries, 0, (size - mt->len_entries) * sizeof(MashEntry));
        
        if ((ptr = realloc(mt->available_entries, size * sizeof(uint32_t))) == NULL) {
            return -1;
            }
        mt->available_entries = (uint32_t *)ptr;
        for (i = mt->len_entries; i < size; ++i) {
            mt->available_entries[i] = i;
            }

        if ((ptr = realloc(mt->keys, size * mt->max_key_len)) == NULL) {
            return -1;
            }
        mt->keys = ptr;

        if ((ptr = realloc(mt->data, size * mt->max_data_len)) == NULL) {
            return -1;
            }
        mt->data = ptr;
        
        mt->len_entries = size;
        }
    else if (size < mt->len_entries && size >= mt->entries_occupied) {
        // TODO  Is this even worth doing?
        }
        
    return 0;
    }



static int resize_keys(MashTable *mt, size_t max_key_len) {
    uint32_t i = 0;
    void *ptr = NULL;
    
    if (max_key_len > mt->max_key_len) {
        if ((ptr = malloc(mt->len_entries * max_key_len)) == NULL) {
            return -1;
            }
        
        for (i = 0; i < mt->len_entries; ++i) { // TODO can be made more efficient - copy only used keys
            memcpy(ptr + (i * max_key_len), mt->keys + (i * mt->max_key_len), mt->max_key_len);
            }
        free(mt->keys);
        mt->keys = ptr;
        mt->max_key_len = max_key_len;
        }
    return 0;
    }



static int resize_data(MashTable *mt, size_t max_data_len) {
    uint32_t i = 0;
    void *ptr = NULL;
    
    if (max_data_len > mt->max_data_len) {
        if ((ptr = malloc(mt->len_entries * max_data_len)) == NULL) {
            return -1;
            }
        
        for (i = 0; i < mt->len_entries; ++i) { // TODO can be made more efficient - copy only used keys
            memcpy(ptr + (i * max_data_len), mt->data + (i * mt->max_data_len), mt->max_data_len);
            }
        free(mt->data);
        mt->data = ptr;
        mt->max_data_len = max_data_len;
        }
    return 0;
    }

