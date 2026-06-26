/*
 * mandelbrot.c - Mandelbrot set generator
 * Copyright (C) 2026 asdftowel
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Based on C gcc #6 by Kevin Miller */

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#if !defined(__STDC_VERSION__) || \
  __STDC_VERSION__ < 201112L || \
  defined(__STDC_NO_THREADS__)
#error "This program requires support for C11 threads."
#elif !defined(UINT8_MAX)
#error "This program requires support for C99 exact width integers."
#else
#include <threads.h>
#endif

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

static inline void calc_sum(
    double * restrict const sums,
    double * restrict const reals,
    double * restrict const imags,
    double const * restrict const init_r,
    double const init_i
) {
  int i;
  double re, im, r2, i2;
  for (i = 0; i < 8; ++i) {
    re = reals[i];
    im = imags[i];
    r2 = re * re;
    i2 = im * im;
    sums[i] = r2 + i2;
    reals[i] = r2 - i2 + init_r[i];
    imags[i] = re * im * 2 + init_i;
  }
}

#define VEC_ALL_GT(vec, cmp) \
  (vec[0] > cmp &&	     \
   vec[1] > cmp &&	     \
   vec[2] > cmp &&	     \
   vec[3] > cmp &&	     \
   vec[4] > cmp &&	     \
   vec[5] > cmp &&	     \
   vec[6] > cmp &&	     \
   vec[7] > cmp)

static inline unsigned int mand8(
    double const * restrict const init_r,
    double const init_i
) {
  double
    reals[8] = {
    init_r[0],
    init_r[1],
    init_r[2],
    init_r[3],
    init_r[4],
    init_r[5],
    init_r[6],
    init_r[7]
  },
    imags[8] = {
      init_i, init_i, init_i, init_i, init_i, init_i, init_i, init_i
    },
    sums[8];
  unsigned int result = 0xff, checks[8];
  int i, j;
  for (i = 0; i < 6; ++i) {
    for (j = 0; j < 8; ++j) {
      calc_sum(sums, reals, imags, init_r, init_i);
    }
    if (VEC_ALL_GT(sums, 4.)) {
      result = 0;
      goto quick_ret;
    }
  }
  calc_sum(sums, reals, imags, init_r, init_i);
  calc_sum(sums, reals, imags, init_r, init_i);
  checks[0] = sums[0] <= 4.;
  checks[1] = sums[1] <= 4.;
  checks[2] = sums[2] <= 4.;
  checks[3] = sums[3] <= 4.;
  checks[4] = sums[4] <= 4.;
  checks[5] = sums[5] <= 4.;
  checks[6] = sums[6] <= 4.;
  checks[7] = sums[7] <= 4.;
  result &= 0x7fu | (checks[0] << 7u);
  result &= 0xbfu | (checks[1] << 6u);
  result &= 0xdfu | (checks[2] << 5u);
  result &= 0xefu | (checks[3] << 4u);
  result &= 0xf7u | (checks[4] << 3u);
  result &= 0xfbu | (checks[5] << 2u);
  result &= 0xfdu | (checks[6] << 1u);
  result &= 0xfeu | checks[7];
 quick_ret:
  return result;
}

#undef VEC_ALL_GT

struct thread_args {
  uint8_t * bitmap;
  double const * init_r;
  double const * init_i;
  unsigned long width;
  unsigned long nrows;
};

static inline int compute8(void *args) {
  struct thread_args const * const locals = args;
  uint8_t * const bitmap = locals->bitmap;
  double const * const init_r = locals->init_r,
    * const init_i = locals->init_i;
  unsigned long const width = locals->width,
    nrows = locals->nrows;
  unsigned long row, col, ir_offset;
  size_t row_offset;
  for (row = 0; row < nrows; ++row) {
    row_offset = width * row;
    for (col = 0, ir_offset = 0; col < width; ++col, ir_offset += 8) {
      bitmap[row_offset + col] = (uint8_t)mand8(init_r + ir_offset, init_i[row]);
    }
  }
  free(args);
  return 0;
}

#ifndef MAX_IMG_SIZE
#define MAX_IMG_SIZE 16000ul /* In pixels */
#endif

#ifndef N_THREADS
#define N_THREADS 4 /* Adjust to number of hardware threads */
#elif N_THREADS < 1
#error "N_THREADS must be at least 1."
#endif

#define SET_THREAD_ARGS(name, off, bm, r, i, w, rows)		\
  name->bitmap = bm + off * w;					\
  name->init_r = r;						\
  name->init_i = i + off;					\
  name->width = w;					       	\
  name->nrows = rows

int main(int argc, char *argv[]) {
  unsigned long img_size, row_size, i, remaining, eq_parts;
  size_t bitmap_size, t, offset;
  uint8_t *bitmap;
  char *end = NULL;
  double inv, *reals, *imags;
  thrd_t threads[N_THREADS];
  struct thread_args * args;
  _Bool failed = 0;

#ifdef _WIN32
  /* Windows silently converts \n into \r\n */
  if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
    perror("Error");
    return EXIT_FAILURE;
  }
#endif
  
  if (argc != 2) {
    fputs("Usage: ./mandelbrot <image size>\n", stderr);
    return EXIT_FAILURE;
  } else if (!(img_size = strtoul(argv[1], &end, 10))) {
    fprintf(stderr, "Error: Zero pixels. Is \"%s\" a number?\n", argv[1]);
    return EXIT_FAILURE;
  } else if (errno == ERANGE) {
    perror("Error");
    return EXIT_FAILURE;
  }

  if (ULONG_MAX - 7ul < img_size) {
    fputs("Warn: Rounding up would cause overflow, rounding down.\n", stderr);
  } else {
    img_size += 7ul;
  }
  
  img_size &= ~7ul;
  inv = 2. / (double)img_size;
  row_size = img_size / 8;
  eq_parts = img_size / N_THREADS;
  remaining = img_size % N_THREADS;
  if (!eq_parts) {
    /*
     * If there's not enough work for all threads,
     * it can probably be done by two.
     */
    eq_parts = img_size / 2;
    remaining = img_size % 2;
  }
  
  if (img_size > MAX_IMG_SIZE) {
    fputs("Error: Image size larger than compile-time limit.", stderr);
    return EXIT_FAILURE;
  }
  
  if (SIZE_MAX / row_size < img_size) {
    fputs(
	"Error: Input dimensions would cause overflow, "
	"try decreasing the size.\n",
	stderr
    );
    return EXIT_FAILURE;
  } else {
    bitmap_size = (size_t)img_size * row_size;
  }

  if (!(bitmap = malloc(bitmap_size))) {
    fputs("Error: cannot allocate bitmap buffer (out of memory?)\n", stderr);
    return EXIT_FAILURE;
  } else if (!(reals = malloc(img_size * sizeof(double)))) {
    fputs("Error: cannot allocate real numbers (out of memory?)\n", stderr);
    free(bitmap);
    return EXIT_FAILURE;
  } else if (!(imags = malloc(img_size * sizeof(double)))) {
    fputs("Error: cannot allocate imaginary numbers (out of memory?)\n", stderr);
    free(bitmap);
    free(reals);
    return EXIT_FAILURE;
  }
  for (i = 0; i < img_size; ++i) {
    reals[i] = inv * (double)i - 1.5;
    imags[i] = inv * (double)i - 1.;
  }
  t = 0;
  offset = 0;
  if (remaining) {
    for (; remaining > 0; ++t, --remaining) {
      if (!(args = malloc(sizeof * args))) {
	goto thread_fail;
      }
      SET_THREAD_ARGS(
	  args,
	  offset,
	  bitmap,
	  reals,
	  imags,
	  row_size,
	  eq_parts + 1
      );
      offset += eq_parts + 1;
      if (thrd_create(threads + t, compute8, args) != thrd_success) {
	goto thread_fail;
      }
    }
  }
  for (; t < N_THREADS; ++t) {
    if (!(args = malloc(sizeof * args))) {
      goto thread_fail;
    }
    SET_THREAD_ARGS(
	args,
	offset,
	bitmap,
	reals,
	imags,
	row_size,
	eq_parts
    );
    offset += eq_parts;
    if (thrd_create(threads + t, compute8, args) != thrd_success) {
      goto thread_fail;
    }
  }
  for (t = 0; t < N_THREADS; ++t) {
    if (thrd_join(threads[t], NULL) != thrd_success) {
      fprintf(stderr, "Error: thread #%zu failed, exiting.\n", t + 1);
      failed = 1;
    }
  }
  free(imags);
  free(reals);
  if (failed) {
    goto late_fail;
  } else if (
     printf("P4\n%lu %lu\n", img_size, img_size) < 0 ||
     fwrite(bitmap, sizeof(uint8_t), bitmap_size, stdout) != bitmap_size
  ) {
    fputs("Error: data output failed, exiting.\n", stderr);
    goto late_fail;
  }
  free(bitmap);
  return EXIT_SUCCESS;
 thread_fail:
  /* Uh oh, let previous threads finish and exit */
  fputs("Error: thread creation failed (out of memory?)\n", stderr);
  if (t) {
    --t;
    for (; t < SIZE_MAX; --t) {
      /* Don't care about thread status because nothing will be written */
      thrd_join(threads[t], NULL);
    }
  }
  free(imags);
  free(reals);
 late_fail:
  free(bitmap);
  return EXIT_FAILURE;
}

#undef N_THREADS
#undef MAX_IMG_SIZE
#undef SET_THREAD_ARGS
