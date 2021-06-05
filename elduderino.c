#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "elduderino.h"
#include "hash.h"
#include "mash.h"


#define NONE 0
#define THRUPLEX 1
#define THRUPLEX_HV 2
#define PRISM 3

const uint16_t UNMAPPED = 0x4;
const uint16_t MATE_UNMAPPED = 0x8;
const uint16_t REVERSE = 0x10;
const uint16_t MATE_REVERSE = 0x20;
const uint16_t READ1 = 0x40;
const uint16_t READ2 = 0x80;
const uint16_t SECONDARY = 0x100;
const uint16_t FILTERED = 0x200;
const uint16_t SUPPPLEMENTARY = 0x800;
const uint16_t NON_PRIMARY = 0x100 | 0x800; // SECONDARY | SUPPPLEMENTARY;
const uint16_t BOTH_UNMAPPED = 0x4 | 0x8; // UNMAPPED | MATE_UNMAPPED;
const uint16_t READXREVERSEXUNMAPPEDX = 0x40 | 0x80 | 0x10 | 0x20 | 0x4 | 0x8; // READ1 | READ2 | REVERSE | MATE_REVERSE | UNMAPPED | MATE_UNMAPPED;
const uint16_t READX = 0x40 | 0x80; // READ1 | READ2;

const char *CONSUMES_REF = "MDN=X";
const char *CONSUMES_READ = "MIS=X";


static bool endswith(const char *text, const char *suffix);
static Segment parse_segment(char *sam, char *sam_end);
static void segment_fprintf(Segment segment, FILE *fp);
static int32_t cigar_len(char *cigar, size_t cigar_len, const char *ops);



int main (int argc, char **argv) {
    /*
     * 
     */
    char *input_filename = NULL, *output_filename = "deduplicated.sam", *stats_filename = "stats.json";
    int min_family_size = 1, umi = NONE;
    
    int sam_fd = -1;
    size_t sam_len = 0;
    char *sam_start = NULL, *sam_end = NULL, *sam = NULL;
    
    char *current_rname = "", *sort_check_rname = "", *mate = NULL, *position = NULL, *next = NULL, *read = NULL;
    size_t current_rname_len = 0, sort_check_rname_len = 0, len = 0, max_position_len = 0, position_len = 0;
    int32_t max_pos = 0, max_pos2 = 0, sort_check_pos = 0, segment_begin = 0, mate_begin = 0;
    HashTable *unpaired = NULL;
    MashTable *paired = NULL, *paired2 = NULL, *swap_mash = NULL;
    Segment segment, mate_segment, swap_segment;
    ReadPair readpair;
    
    char *data = NULL, *key = NULL;
    size_t data_size = 0, key_size = 0;
    uint32_t bucket = 0;
    
    
    
    char *endptr = NULL;
    long val = 0;
    // Variables needed by getopt_long
    int option_index = 0, c = 0;
    static struct option long_options[] = {{"output", required_argument, 0, 'o'},
                                           {"stats", required_argument, 0, 's'},
                                           {"umi", required_argument, 0, 'u'},
                                           {"min-family-size", no_argument, 0, 'm'},
                                           {0, 0, 0, 0}};

    // Parse optional arguments
    while (c != -1) {
        c = getopt_long (argc, argv, "o:s:e:npa:m:l:u:", long_options, &option_index);

        switch (c) {
            case 'o':
                if (!endswith(optarg, ".sam")) {
                    fprintf(stderr, "Error: Output file must be of type sam\n");
                    exit(EXIT_FAILURE);
                    }
                output_filename = optarg;
                break;
                
            case 's':
                if (!endswith(optarg, ".json")) {
                    fprintf(stderr, "Error: Stats file must be of type json\n");
                    exit(EXIT_FAILURE);
                    }
                stats_filename = optarg;
                break;
            
            case 'm':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg)) {
                    fprintf(stderr, "Error: Invalid --min-family-size\n");
                    exit(EXIT_FAILURE);
                    }
                min_family_size = (int)val;
                break;                
                
            case 'u':
                if (strcmp(optarg, "thruplex")) {
                    umi = THRUPLEX;
                    }
               else if (strcmp(optarg, "thruplex_hv")) {
                    umi = THRUPLEX_HV;
                    }
                else if (strcmp(optarg, "prism")) {
                    umi = PRISM;
                    }
                else {
                    fprintf(stderr, "Error: Unsupported umi type: %s\n", optarg);
                    exit(EXIT_FAILURE);
                    }
                break;
                
            case '?':
                // unknown option, getopt_long already printed an error message.
                exit(EXIT_FAILURE);
            }
        }
    
    
    if (argc - optind != 1) {
        fprintf(stderr, "Error: No input file supplied.\n");
        exit(EXIT_FAILURE);
        }
    input_filename = *(argv + optind);
    if (!endswith(input_filename, ".sam")) {
        fprintf(stderr, "Error: Input file must be of type sam\n");
        exit(EXIT_FAILURE);
        }
    
    
    if ((sam_fd = open(input_filename, O_RDONLY)) == -1) {
        fprintf(stderr, "Error: Unable to open %s\n", input_filename);
        exit(EXIT_FAILURE);
        }
    
    sam_len = (size_t)lseek(sam_fd, 0, SEEK_END);
    lseek(sam_fd, 0, SEEK_SET);
    if ((sam_start = mmap(NULL, sam_len, PROT_READ, MAP_PRIVATE, sam_fd, 0)) == NULL) {
        fprintf(stderr, "Error: Unable to memory map sam file\n");
        exit(EXIT_FAILURE);
        }
    
    
    // Move past all comments to first read
    sam_end = sam_start + sam_len;
    for (sam = sam_start; sam < sam_end; ++sam) {
        if (*sam != '@') {
            break;
            }
        for (; *sam != '\n'; ++sam);
        }
    
    unpaired = hash_new(64);
    paired = mash_new(64);
    paired2 = mash_new(64);
        
    for (;sam < sam_end; sam = next) {
        segment = parse_segment(sam, sam_end);
        next = segment.next;
        // Sanity check to ensure that sam file is sorted by position
        if (segment.rname_len == sort_check_rname_len && memcmp(segment.rname, sort_check_rname, sort_check_rname_len) == 0) {
            if (segment.pos < sort_check_pos) {
                fprintf(stderr, "Error: Sam file must be sorted by position\n");
                exit(EXIT_FAILURE);
                }
            }
        else {
            sort_check_rname = segment.rname;
            sort_check_rname_len = segment.rname_len;
            }
        sort_check_pos = segment.pos;
        
        // Skip secondary, supplementary and completely unmapped reads
        if ((segment.flag & NON_PRIMARY) || ((segment.flag & BOTH_UNMAPPED) == BOTH_UNMAPPED)) {
            continue;
            }
        
        //-------------------------------------------hash_validate(unpaired, stderr);
        
        // Do we have a pair of reads yet? If not store this read and move on to the next
        if ((mate = hash_pop(unpaired, segment.qname, segment.qname_len, &len)) == NULL) {
            //------------------------------------------hash_summary_fprintf(unpaired, stderr);
            hash_put(unpaired, segment.qname, segment.qname_len, segment.qname, segment.len);
            continue;
            }
        
        mate_segment = parse_segment(mate, mate + len);
        
        // This is needed as an unmapped read may be positioned before or after its mate depending
        // on the value of the REVERSE flag in a sorted sam
        // If one read is unmapped then both reads will share the same rname and pos therefore these
        // don't get affected by the swap
        if (mate_segment.flag & UNMAPPED) {
            swap_segment = segment;
            segment = mate_segment;
            mate_segment = swap_segment;
            }
        
        
        // Segment_begin and mate_begin are the positions of the start of a read which will be the
        // rightmost position if reverse complemnted
        mate_begin = mate_segment.pos;
        if (mate_segment.flag & REVERSE) {
            mate_begin += cigar_len(mate_segment.cigar, mate_segment.cigar_len, CONSUMES_REF);
            }
        
        if (segment.flag & UNMAPPED) {
            segment_begin = mate_begin;
            }
        else {
            segment_begin = segment.pos;
            if (segment.flag & REVERSE) {
                segment_begin += cigar_len(segment.cigar, segment.cigar_len, CONSUMES_REF);
                }
            }
        
        
        // 2 x max number of decimal digits in an int32 (10 digits) + 
        // max number of decimal digits in a uint16 (5 digits) +
        // 4 x '\t' field separators +
        // terminal '\0'
        position_len = mate_segment.rname_len + segment.rname_len + 30;
        if (position_len > max_position_len) {
            if ((position = realloc(position, position_len)) == NULL) {
                fprintf(stderr, "Error: Unable to alllocate memory for position buffer\n");
                exit(EXIT_FAILURE);
                }
            max_position_len = position_len;
            }
        snprintf(position, position_len, "%.*s\t%010i\t%.*s\t%010i\t%05u",
                                         (int)mate_segment.rname_len, mate_segment.rname,
                                         (int)mate_begin,
                                         (int)segment.rname_len, segment.rname,
                                         (int)segment_begin,
                                         (unsigned)(segment.flag & READXREVERSEXUNMAPPEDX));
        
        readpair.left = mate_segment;
        readpair.right = segment;
        
        if (mate_begin > segment_begin) {
            segment_begin = mate_begin;
            }
        
        if (segment.rname_len == current_rname_len && memcmp(segment.rname, current_rname, current_rname_len) == 0 && segment.pos <= max_pos) {
            if (mash_get(paired, position, position_len, &len) != NULL) {
                if (mash_put(paired, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                    fprintf(stderr, "Error: Unable to to add position to paired hash table\n");
                    exit(EXIT_FAILURE);
                    }
                }
            else {
                if (mash_put(paired2, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                    fprintf(stderr, "Error: Unable to to add position to paired2 hash table\n");
                    exit(EXIT_FAILURE);
                    }
                if (segment_begin > max_pos2) {
                    max_pos2 = segment_begin;
                    }
                }
            }
        else {
            bucket = 0;
            key = NULL;
            while ((data = mash_popall(paired, (const void **)&key, &key_size, &data_size, &bucket)) != NULL) {
                //fprintf(stderr, "%.*s\n", (int)key_size, key);//, (int)data_size, data);
                }
            swap_mash = paired;
            paired = paired2;
            max_pos = max_pos2;
            paired2 = swap_mash;
            max_pos2 = 0;
            if (mash_put(paired, position, position_len, &readpair, sizeof(ReadPair)) == -1) {
                fprintf(stderr, "Error: Unable to to add position to paired hash table\n");
                exit(EXIT_FAILURE);
                }
            if (segment.rname_len != current_rname_len || memcmp(segment.rname, current_rname, current_rname_len) != 0) {
                current_rname = segment.rname;
                current_rname_len = segment.rname_len;
                max_pos = segment_begin;
                }
            else if (segment_begin > max_pos) {
                max_pos = segment_begin;
                }
            }
        
        }
    
    //mash_summary_fprintf(paired, stderr);
    //mash_summary_fprintf(paired2, stderr);
    }



void dedupe(ReadPair *family, size_t family_size, int min_family_size) {
    }



static Segment parse_segment(char *read, char *sam_end) {
    Segment segment;
    int column = 0;
    char *start = NULL, *endptr = NULL, *beginning = read;
    long val = 0;
    
    segment.start = read;
    segment.barcode = NULL; // optional field so must be initialised
    segment.barcode_len = 0; // optional field so must be initialised
    for (start = read; read < sam_end; ++read) {
        if (*read == '\t' || *read == '\n') {
            switch (++column) {
                case 1: // qname
                    segment.qname = start;
                    segment.qname_len = read - start;
                    break;
                case 2: // flag
                    errno = 0;
                    val = strtol(start, &endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid flag in sam file.\n");
                        exit(EXIT_FAILURE);
                        }
                    segment.flag = (uint16_t)val;
                    break;
                case 3: // rname
                    segment.rname = start;
                    segment.rname_len = read - start;
                    break;
                case 4: // pos
                    errno = 0;
                    val = strtol(start, &endptr, 10);
                    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == start)) {
                        fprintf(stderr, "Error: Invalid pos in sam file.\n");
                        exit(EXIT_FAILURE);
                        }
                    segment.pos = (int32_t)val;
                    break;
                case 6: // cigar
                    segment.cigar = start;
                    segment.cigar_len = read - start;
                    break;
                case 10: // seq
                    segment.seq = start;
                    segment.seq_len = read - start;
                    break;
                case 11: // qual
                    segment.qual = start;
                    segment.qual_len = read - start;
                    break;
                default:
                    if (column > 11) { // barcode tag
                        if (memcmp(start, "RX:Z:", 5) == 0) {
                            segment.barcode = start + 5;
                            segment.barcode_len = read - start - 5;
                            }
                        }
                }
            start = read + 1;
            if (*read == '\n') {
                segment.len = read - beginning + 1; // includes final \n
                segment.next = read + 1;
                break;
                }
            }
        }
    if (column < 11) {
        fprintf(stderr, "Error: Truncated sam file.\n");
        exit(EXIT_FAILURE);
        }
    return segment;
    }



static bool endswith(const char *text, const char *suffix) {
    int offset = strlen(text) - strlen(suffix);
    return (offset >= 0 && strcmp(text + offset, suffix) == 0);
    }



static int32_t cigar_len(char *cigar, size_t cigar_len, const char *ops) {
    long val = 0;
    int32_t len = 0;
    char *endptr = NULL, *cigar_end = cigar + cigar_len;

    if (*cigar == '*') {
        return 0;
        }
    
    for (; cigar < cigar_end; cigar = endptr + 1) {
        errno = 0;
        val = strtol(cigar, &endptr, 10);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == cigar)) {
            fprintf(stderr, "Error: Invalid cigar string %.*s\n", (int)cigar_len, cigar);
            exit(EXIT_FAILURE);
            }
            
        if (strchr(ops, *endptr) != NULL) {
            len += (int32_t)val;
            }
        }
    return len;
    }



static void segment_fprintf(Segment segment, FILE *fp) {
    fprintf(fp, "QNAME: %.*s\n", (int)segment.qname_len, segment.qname);
    fprintf(fp, "FLAG: ");
    if (segment.flag & UNMAPPED) fprintf(fp, " UNMAPPED");
    if (segment.flag & MATE_UNMAPPED) fprintf(fp, " MATE_UNMAPPED");
    if (segment.flag & REVERSE) fprintf(fp, " REVERSE");
    if (segment.flag & MATE_REVERSE) fprintf(fp, " MATE_REVERSE");
    if (segment.flag & READ1) fprintf(fp, " READ1");
    if (segment.flag & READ2) fprintf(fp, " READ2");
    if (segment.flag & SECONDARY) fprintf(fp, " SECONDARY");
    if (segment.flag & FILTERED) fprintf(fp, " FILTERED");
    if (segment.flag & SUPPPLEMENTARY) fprintf(fp, " SUPPPLEMENTARY");
    fprintf(fp, "\n");
    fprintf(fp, "RNAME: %.*s\n", (int)segment.rname_len, segment.rname);
    fprintf(fp, "POS:   %i\n", (int)segment.pos);
    fprintf(fp, "CIGAR: %.*s\n", (int)segment.cigar_len, segment.cigar);
    }

