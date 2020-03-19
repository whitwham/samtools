/*  bam_clip.c -- clip reads.

    Copyright (C) 2020 Genome Research Ltd.

    Author: Andrew Whitwham <aw7@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE
*/

#include <config.h>

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include <htslib/hts.h>
#include "htslib/hfile.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "samtools.h"

typedef struct {
    int64_t left;
    int64_t right;
} bed_pair_t;



static void print_bed_pairs(bed_pair_t *bp, int len) {
    int i;
    
    for (i = 0; i < len; i++) {
        printf("%ld %ld\n", bp[i].left, bp[i].right);
    }
}


static bed_pair_t *load_bed_file_pairs(char *infile, int *length) {
    hFILE *fp;
    bed_pair_t *bp;
    int size = 256, len = 0, line_count = 0;
    int64_t left, right;
    kstring_t line = {0, 0, NULL};
    
    if ((fp = hopen(infile, "r")) == NULL) {
        fprintf(stderr, "[clip] error: unable to open file %s\n", infile);
        return NULL;
    }
    
    if ((bp = malloc(size * sizeof(bed_pair_t))) == NULL) {
        fprintf(stderr, "[clip] error: unable to allocate memory for bed data.\n");
        goto error;
    }
    
    while (line.l = 0, kgetline(&line, (kgets_func *)hgets, fp) >= 0) {
        line_count++;
    
        if (line.l == 0 || *line.s == '#') continue;
        if (strncmp(line.s, "track ", 6) == 0) continue;
        if (strncmp(line.s, "browser ", 8) == 0) continue;
        
        if (sscanf(line.s, "%*s %ld %ld", &left, &right) != 2) {
            fprintf(stderr, "[clip] error: bad bed file format in line %d of %s",
                                line_count, infile);
            free(bp);
            bp = NULL;
            goto error;
        }
        
        if (len == size) {
            bed_pair_t *tmp;
           
            size *= 2;
           
            if ((tmp = realloc(bp, size * sizeof(bed_pair_t))) == NULL) {
                fprintf(stderr, "[clip] error: unable to allocate more memory for bed data.\n");
                free(bp);
                bp = NULL;
                goto error;
            }
            
            bp = tmp;
        }
        
        bp[len].left  = left;
        bp[len].right = right;
        len++;
    }
    
    *length = len;
    
error:        
    ks_free(&line);
    if (hclose(fp) != 0) {
        fprintf(stderr, "[clip] Warning: failed to close %s", infile);
    }
    
    return bp;
}
    
    
static int matching_clip_site(bed_pair_t *sites, int s_length, hts_pos_t pos) {
    int i, tol = 5, size;  // may need this to be variable
    
    // fprintf(stderr, "[debug] called for s_length %d pos %ld\n", s_length, pos);
    
    for (i = 0; i < s_length; i++) {
        hts_pos_t mod_left, mod_right;
    
        size = 0;
        
        if (sites[i].left > tol) {
            mod_left = sites[i].left - tol;
        } else {
            mod_left = 0;
        }
        
        mod_right = sites[i].right + tol;
        
        if (pos >= mod_left && pos <= mod_right) {
            size = sites[i].right - sites[i].left;
            break;
        }
    }
    
    return size;
}


static int bam_trim_left(bam1_t *rec, uint32_t bases) {
    uint32_t *orig_cigar = bam_get_cigar(rec);
    uint8_t *orig_seq = bam_get_seq(rec);
    uint8_t *orig_qual = bam_get_qual(rec);
    uint8_t *orig_aux = bam_get_aux(rec);
    uint32_t *new_cigar;
    uint32_t new_n_cigar;
    uint8_t *new_qual;
    size_t orig_l_aux = bam_get_l_aux(rec);
    uint32_t i, j;
    uint32_t removed = bases, hardclip = 0;
    hts_pos_t new_pos = rec->core.pos;
    uint32_t cig_type;
    uint8_t *new_seq_start;
    
    if (bases >= rec->core.l_qseq) {
        rec->core.l_qseq = 0;
        rec->core.n_cigar = 0;
        if (orig_l_aux)
            memmove(bam_get_aux(rec), orig_aux, orig_l_aux);
        rec->l_data -= orig_aux - bam_get_aux(rec);
        return 1;
    }

    // Shrink CIGAR here
    
    // make memory for a temporary cigar string, one bigger than the original
    // to account for the extra H operation
    if ((new_cigar = malloc((rec->core.n_cigar + 1) * sizeof(uint32_t))) == NULL) {
        fprintf(stderr, "[clip] error: unable to allocate memory for new cigar string.\n");
        return 1;
    }
    
    for (i = 0;  i < rec->core.n_cigar; i++) {
        cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));
        
        if (cig_type == 0) {
            hardclip += bam_cigar_oplen(orig_cigar[i]);
        } else {
            if (cig_type & 1) {
                if (bam_cigar_oplen(orig_cigar[i]) < removed) {
                    removed -= bam_cigar_oplen(orig_cigar[i]);
                } else {
                    break;
                }
            }
            
            if (cig_type & 2) {
                new_pos += bam_cigar_oplen(orig_cigar[i]);
            }
        }
    }
    
    cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));
    
    // account for the last operation
    if (cig_type & 2) {
        new_pos += removed;
    }

    j = 0;
    new_cigar[j++] = bam_cigar_gen(hardclip + bases, BAM_CHARD_CLIP);

    if (bam_cigar_oplen(orig_cigar[i]) != removed) {
        new_cigar[j++] = bam_cigar_gen(bam_cigar_oplen(orig_cigar[i]) - removed, bam_cigar_op(orig_cigar[i]));
    }

    // fill in the rest of the cigar
    i++;
    
    for (; i < rec->core.n_cigar; i++) {
         new_cigar[j++] = orig_cigar[i];
    }
    
    new_n_cigar = j;
    
    if (new_n_cigar != rec->core.n_cigar) {
        new_seq_start = orig_seq + ((new_n_cigar - rec->core.n_cigar) * sizeof(uint32_t));
    } else {
        new_seq_start = orig_seq;
    }
    
    new_qual = new_seq_start + (rec->core.l_qseq - bases + 1) / 2;
    // Move / shrink SEQ
    if ((bases & 1) == 0) {
        memmove(new_seq_start, orig_seq + (bases / 2),
	        (rec->core.l_qseq - bases) / 2);
    } else {
        uint8_t *in = orig_seq + bases / 2;
        uint8_t *out = new_seq_start;
        uint32_t i;
        for (i = bases; i < rec->core.l_qseq - 1; i += 2) {
            *out++ = ((in[0] & 0x0f) << 4) | ((in[1] & 0xf0) >> 4);
            in++;
        }
        if (i < rec->core.l_qseq) {
            *out++ = (in[0] & 0x0f) << 4;
        }
        assert(out == new_qual);
    }
    
    // now the sequence is in position put the new cigar string in
    memcpy(orig_cigar, new_cigar, new_n_cigar * sizeof(uint32_t));
    rec->core.n_cigar = new_n_cigar;

    // Move / shrink QUAL
    memmove(new_qual, orig_qual, rec->core.l_qseq - bases);

    // Set new l_qseq
    rec->core.l_qseq -= bases;

    // Move AUX
    if (orig_l_aux)
        memmove(bam_get_aux(rec), orig_aux, orig_l_aux);

    // Set new l_data
    rec->l_data -= orig_aux - bam_get_aux(rec);
    
    // put in new pos
    rec->core.pos = new_pos;
    
    return 0;
}


static int bam_trim_right(bam1_t *rec, uint32_t bases) {
    uint32_t *orig_cigar = bam_get_cigar(rec);
    uint8_t *orig_seq = bam_get_seq(rec);
    uint8_t *orig_qual = bam_get_qual(rec);
    uint8_t *orig_aux = bam_get_aux(rec);
    uint32_t *new_cigar;
    uint32_t new_n_cigar = 0;
    uint8_t *new_qual;
    size_t orig_l_aux = bam_get_l_aux(rec);
    uint32_t i;
    int32_t j;
    uint32_t removed = bases, hardclip = 0;
    uint32_t cig_type;
    uint8_t *new_seq_start;

    if (bases >= rec->core.l_qseq) {
        rec->core.l_qseq = 0;
        rec->core.n_cigar = 0;
        if (orig_l_aux)
            memmove(bam_get_aux(rec), orig_aux, orig_l_aux);
        rec->l_data -= orig_aux - bam_get_aux(rec);
        return 1;
    }
    
    // Shrink CIGAR here

    // make memory for a temporary cigar string, one bigger than the original
    // to account for the extra H operation
    if ((new_cigar = malloc((rec->core.n_cigar + 1) * sizeof(uint32_t))) == NULL) {
        fprintf(stderr, "[clip] error: unable to allocate memory for new cigar string.\n");
        return 1;
    }

    for (i = rec->core.n_cigar - 1;  i > 0; i--) {
        cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));
        
        if (cig_type == 0) {
            hardclip += bam_cigar_oplen(orig_cigar[i]);
        } else {
            if (cig_type & 1) {
                if (bam_cigar_oplen(orig_cigar[i]) < removed) {
                    removed -= bam_cigar_oplen(orig_cigar[i]);
                } else {
                    break;
                }
            }
        }
    }
    
    cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));
    
    j = rec->core.n_cigar;
    new_cigar[j--] = bam_cigar_gen(hardclip + bases, BAM_CHARD_CLIP);
    new_n_cigar++;

    if (bam_cigar_oplen(orig_cigar[i]) != removed) {
        new_cigar[j--] = bam_cigar_gen(bam_cigar_oplen(orig_cigar[i]) - removed, bam_cigar_op(orig_cigar[i]));
        new_n_cigar++;
    }
    
    // fill in the rest of the cigar
    i++;
    
    for (; i < rec->core.n_cigar; i++) {
        new_cigar[j--] = orig_cigar[i];
        new_n_cigar++;
    }
    
    if (new_n_cigar != rec->core.n_cigar) {
        new_seq_start = orig_seq + ((new_n_cigar - rec->core.n_cigar) * sizeof(uint32_t));
    } else {
        new_seq_start = orig_seq;
    }

    new_qual = new_seq_start + (rec->core.l_qseq - bases + 1) / 2;
    // Move SEQ
    memmove(new_seq_start, orig_seq, (rec->core.l_qseq - bases + 1) / 2);

    // now the sequence is in position put the new cigar string in
    memcpy(orig_cigar, new_cigar, new_n_cigar * sizeof(uint32_t));
    rec->core.n_cigar = new_n_cigar;

    // Move / shrink QUAL
    memmove(new_qual, orig_qual, rec->core.l_qseq - bases);

    // Set new l_qseq
    rec->core.l_qseq -= bases;

    // Move AUX
    if (orig_l_aux)
        memmove(bam_get_aux(rec), orig_aux, orig_l_aux);

    // Set new l_data
    rec->l_data -= orig_aux - bam_get_aux(rec);
}


static int bam_clip(samFile *in, samFile *out, char *bedfile, int add_pg, char *args) {
    bed_pair_t *positions;
    int pos_length, ret = 0, exclude = 0;
    bam_hdr_t *header;
    bam1_t *b;
    long f_count = 0, r_count = 0, n_count = 0, l_count = 0, l_exclude = 0;
    
    if ((positions = load_bed_file_pairs(bedfile, &pos_length)) == NULL) {
        ret = 1;
        goto fail;
    }
    
    // print_bed_pairs(positions, pos_length);

    if ((header = sam_hdr_read(in)) == NULL) {
        fprintf(stderr, "[clip] error reading header\n");
        ret = 1;
        goto fail;
    }

    if (add_pg && sam_hdr_add_pg(header, "samtools", "VN", samtools_version(),
                        args ? "CL" : NULL,
                        args ? args : NULL,
                        NULL) != 0) {
        fprintf(stderr, "[clip] warning: unable to add @PG line to header.\n");
    }
    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "[clip] error writing header.\n");
        goto fail;
    }
    
    // TODO add write index
    b = bam_init1();
    
    while ((ret = sam_read1(in, header, b)) >= 0) {
        hts_pos_t pos;
        int is_rev;
        int p_size;
        
        l_count++;
        
        // exclude |= (BAM_FUNMAP | BAM_FQCFAIL);
        
        if (!exclude) {

            if (bam_is_rev(b)) {
                pos = bam_endpos(b);
                is_rev = 1;
            } else {
                pos = b->core.pos;
                is_rev = 0;
            }

            if ((p_size = matching_clip_site(positions, pos_length, pos))) {
                if (is_rev) {
                    bam_trim_right(b, p_size);
                    r_count++;
                } else {
                    bam_trim_left(b, p_size);
                    f_count++;
                }
            } else {
                n_count++;
            }
        } else {
            l_exclude++;
        }
        
	if (sam_write1(out, header, b) < 0) {
    	    fprintf(stderr, "[clip] ERROR: could not write line %ld.\n", l_count);
	    goto fail;
    	}
    }
   
   
    fprintf(stdout, "f_count %ld r_count %ld n_count %ld l_count %ld, l_exclude %ld\n",
                f_count, r_count, n_count, l_count, l_exclude);
    







fail: // better error handling later
    free(positions);
    return ret;
}


static void usage(void) {
    fprintf(stderr, "Usage: samtools clip -b bedfile <input.bam> <output.bam>\n");
}


int clip_main(int argc, char **argv) {
    int c, ret, add_pg = 1;
    char wmode[3] = {'w', 'b', 0};
    char *bedfile = NULL, *arg_list;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};
    samFile *in = NULL, *out = NULL;
    
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1002},
        {NULL, 0, NULL, 0}
    };
    
    while ((c = getopt_long(argc, argv, "b:@:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'b': bedfile = optarg; break;
            case 1002: add_pg = 0; break;
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': usage(); exit(1);
        }
    }
    
    if (!bedfile) {
        usage();
        exit(1);
    }
    
    if (optind + 2 > argc) {
        usage();
        exit(1);
    }
    
    if ((in = sam_open_format(argv[optind], "rb", &ga.in)) == NULL) {
        print_error_errno("clip", "cannot open input file");
        exit(1);
    }
    sam_open_mode(wmode+1, argv[optind+1], NULL);
    
    if ((out = sam_open_format(argv[optind+1], wmode, &ga.out)) == NULL) {
        print_error_errno("clip", "cannot open output file");
        exit(1);;
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "[clip] error: cannot create thread pool\n");
            exit(1);
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }
    
    if (add_pg) {
        arg_list = stringify_argv(argc + 1, argv - 1);
    } else {
        arg_list = NULL;
    }
    
    ret = bam_clip(in, out, bedfile, add_pg, arg_list);
    
    // cleanup
    sam_close(in);
    
    if (sam_close(out) < 0) {
        fprintf(stderr, "[clip] error:  errorwhile closing output file %s\n", argv[optind+1]);
        exit(1);
    }

    if (p.pool) hts_tpool_destroy(p.pool);
    
    sam_global_args_free(&ga);
    free(arg_list);
        
    exit(ret);
}    
    
   
   
        
