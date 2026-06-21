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
#include <stdio.h>

#define N_THREADS 4 /* Adjust to number of hardware threads */
#if N_THREADS < 1
#error "N_THREADS must be at least 1"
#endif

#if !defined(__STDC_VERSION__) || \
  __STDC_VERSION__ < 201112L || \
  defined(__STDC_NO_THREADS__) || \
  N_THREADS == 1
#define NO_THREADS
#else
#include <threads.h>
#endif

#if __STDC_VERSION__ < 202311L
#include <stdbool.h>
#endif

#if __STDC_VERSION__ >= 202311L
#define UNREACHABLE unreachable()
#elif defined(_MSC_VER)
#define UNREACHABLE __assume(0)
#elif defined(__GNUC__)
#define UNREACHABLE __builtin_unreachable()
#else
#define UNREACHABLE abort()
#endif

#define STRINGIFY(x) #x
#define EXPAND(x) STRINGIFY(x)

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
  for (i = 0; i < CHAR_BIT; ++i) {
    re = reals[i];
    im = imags[i];
    r2 = re * re;
    i2 = im * im;
    sums[i] = r2 + i2;
    reals[i] = r2 - i2 + init_r[i];
    imags[i] = re * im * 2 + init_i;
  }
}

static inline bool all_gt(double const * restrict const vec) {
  int i;
  for (i = 0; i < CHAR_BIT; ++i) {
    if (vec[i] <= 4.) {
      return false;
    }
  }
  return true;
}

enum bit_consts {
  ULONG_BITS = CHAR_BIT * sizeof(unsigned long),
  HIGH_SHIFT = ULONG_BITS - CHAR_BIT,
  TOP_BIT = CHAR_BIT - 1
};

static inline unsigned int mand_char(
    double const * restrict const init_r,
    double const init_i,
    unsigned int const * restrict const masks,
    bool const prune
) {
  double reals[CHAR_BIT], imags[CHAR_BIT], sums[CHAR_BIT];
  unsigned int result = UCHAR_MAX;
  int i, j;

  for (i = 0; i < CHAR_BIT; ++i) {
    reals[i] = init_r[i];
    imags[i] = init_i;
  }

  if (prune) {
    for (i = 0; i < 6; ++i) {
      for (j = 0; j < 8; ++j) {
	calc_sum(sums, reals, imags, init_r, init_i);
      }
      if (all_gt(sums)) {
	result = 0;
	goto done;
      }
    }
    calc_sum(sums, reals, imags, init_r, init_i);
    calc_sum(sums, reals, imags, init_r, init_i);
  } else {
    for (i = 0; i < 50; ++i) {
      calc_sum(sums, reals, imags, init_r, init_i);
    }
  }
  for (i = 0, j = TOP_BIT; i < CHAR_BIT; ++i, --j) {
    result &= masks[i] | ((unsigned int)(sums[i] <= 4.) << j);
  }
 done:
  return result;
}

static inline unsigned long mand_long(
    double const * restrict const init_r,
    double const init_i,
    unsigned int const * restrict const masks,
    bool const prune
) {
  unsigned long result = 0;
  int offset = 0;
  size_t i;
  for (i = 0; i < sizeof(unsigned long); ++i) {
    result =
      (result >> CHAR_BIT) |
      (mand_char(init_r + offset, init_i, masks, prune) << HIGH_SHIFT);
    offset += CHAR_BIT;
  }
  return result;
}

enum batch_size { BYTE, WORD };

union bitmap_t { unsigned char * bytes; unsigned long * words; };

struct thread_args {
  union bitmap_t bitmap;
  enum batch_size batch;
  double const * init_r;
  double const * init_i;
  unsigned long width;
  unsigned long nrows;
  unsigned long top_mask;
  unsigned int const * masks;
};

static int compute(void *args) {
  struct thread_args const * const locals = args;
  union bitmap_t const bitmap = locals->bitmap;
  double const * const init_r = locals->init_r,
    * const init_i = locals->init_i;
  unsigned long const width = locals->width,
    nrows = locals->nrows,
    top_mask = locals->top_mask;
  unsigned int const * const masks = locals->masks;
  unsigned long row, col, ir_offset;
  size_t row_offset = 0;
  bool prune = false;
  switch (locals->batch) {
  case BYTE:
    for (row = 0; row < nrows; ++row) {
      for (col = 0, ir_offset = 0; col < width; ++col, ir_offset += CHAR_BIT) {
	prune = (
	    bitmap.bytes[row_offset + col] =
	    (unsigned char)mand_char(init_r + ir_offset, init_i[row], masks, prune)
	) == 0;
      }
      row_offset += width;
    }
    break;
  case WORD:
    for (row = 0; row < nrows; ++row) {
      for (col = 0, ir_offset = 0; col < width; ++col, ir_offset += ULONG_BITS) {
	prune = ((
	    bitmap.words[row_offset + col] =
	    mand_long(init_r + ir_offset, init_i[row], masks, prune)
	) & top_mask) == 0;
      }
      row_offset += width;
    }
    break;
  default:
    UNREACHABLE;
  }
  free(args);
  return 0;
}

#define MAX_IMG_SIZE 16000ul /* In pixels */
#define SET_THREAD_ARGS(name, off, bm, s, r, i, w, rows, tm, m)	\
  switch (s) {							\
  case BYTE:							\
    name->bitmap.bytes = bm.bytes + off * w;			\
    break;							\
  case WORD:							\
    name->bitmap.words = bm.words + off * w;			\
    break;							\
  default:							\
    UNREACHABLE;						\
  }								\
  name->batch = s;						\
  name->init_r = r;						\
  name->init_i = i + off;					\
  name->width = w;					       	\
  name->nrows = rows;						\
  name->top_mask = tm;						\
  name->masks = m
#define PUT_ERR(msg)						\
  fputs(__FILE__ ", line " EXPAND(__LINE__) ": Error: " msg "\n", stderr)
#define PUT_WARN(msg)						\
  fputs(__FILE__ ", line " EXPAND(__LINE__) ": Warning: " msg "\n", stderr)
#define PRINTF_ERR(fmt, ...)						\
  fprintf(								\
      stderr, 								\
      __FILE__ ", line " EXPAND(__LINE__) ": Error: " fmt "\n",		\
      __VA_ARGS__							\
  )
#define PUT_ERRNO perror(__FILE__ ", line " EXPAND(__LINE__) ": Error")

/*
 * In a real program the NO_THREADS ifdefs should instead be functions
 * in separate files that get conditionally included based on the macro.
 */

int main(int argc, char *argv[]) {
  unsigned long const top_mask = (unsigned long)UCHAR_MAX << HIGH_SHIFT;
  unsigned long img_size, img_rem, row_size, i;
  unsigned int masks[CHAR_BIT];
  int idx, shift;
  size_t bitmap_size;
  union bitmap_t bitmap;
  char * end = NULL;
  double inv, * reals, * imags;
#ifndef NO_THREADS
  unsigned long eq_parts, remaining;
  size_t t, offset;
  thrd_t threads[N_THREADS];
  bool failed = false;
#endif
  struct thread_args * args;
  bool batch_long;
  enum batch_size batch;

  for (idx = 0, shift = TOP_BIT; idx < CHAR_BIT; ++idx, --shift) {
    masks[idx] = ~(1u << shift);
  }

#ifdef _WIN32
  /* Windows silently converts \n into \r\n */
  if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
    PUT_ERRNO;
    return EXIT_FAILURE;
  }
#endif
  
  if (argc != 2) {
    fputs("Usage: ./mandelbrot <image size>\n", stderr);
    return EXIT_FAILURE;
  } else if (!(img_size = strtoul(argv[1], &end, 10))) {
    PRINTF_ERR("Zero pixels. Is \"%s\" a number?", argv[1]);
    return EXIT_FAILURE;
  } else if (errno == ERANGE) {
    PUT_ERRNO;
    return EXIT_FAILURE;
  }

  if ((img_rem = img_size % CHAR_BIT)) {
    img_size -= img_rem;
    if (ULONG_MAX - CHAR_BIT >= img_size) {
      img_size += CHAR_BIT;
    } else {
      PUT_WARN("Rounding up would cause overflow, rounding down.");
    }
  }
  
  batch_long = sizeof(unsigned char) != sizeof(unsigned long) &&
    img_size % ULONG_BITS == 0;
  inv = 2. / img_size;
  row_size = img_size / CHAR_BIT;
#ifndef NO_THREADS
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
#endif
  if (img_size > MAX_IMG_SIZE) {
    PUT_ERR("Image size larger than compile-time limit.");
    return EXIT_FAILURE;
  }
  
  if (SIZE_MAX / row_size < img_size) {
    PUT_ERR("Input dimensions would cause overflow, try decreasing the size.");
    return EXIT_FAILURE;
  } else {
    bitmap_size = (size_t)img_size * row_size;
  }

  if (SIZE_MAX / sizeof(double) < img_size) {
    PUT_ERR("Allocation sizes would cause overflow, try decreasing the size.");
    return EXIT_FAILURE;
  } else if (!(reals = malloc(img_size * sizeof(double)))) {
    PUT_ERR("Cannot allocate real numbers (out of memory?)");
    return EXIT_FAILURE;
  } else if (!(imags = malloc(img_size * sizeof(double)))) {
    PUT_ERR("Cannot allocate imaginary numbers (out of memory?)");
    free(reals);
    return EXIT_FAILURE;
  }

  for (i = 0; i < img_size; ++i) {
    reals[i] = inv * i - 1.5;
    imags[i] = inv * i - 1.;
  }

  if (batch_long) {
    batch = WORD;
    row_size /= sizeof(unsigned long);
    if (!(bitmap.words = malloc(bitmap_size))) {
      PUT_ERR("Cannot allocate bitmap (out of memory?)");
      free(imags);
      free(reals);
      return EXIT_FAILURE;
    }
    bitmap_size /= sizeof(unsigned long);
  } else {
    batch = BYTE;
    if (!(bitmap.bytes = malloc(bitmap_size))) {
      PUT_ERR("Cannot allocate bitmap (out of memory?)");
      free(imags);
      free(reals);
      return EXIT_FAILURE;
    }
  }
#ifdef NO_THREADS
  if (!(args = malloc(sizeof * args))) {
    free(imags);
    free(reals);
    batch_long ? free(bitmap.words) : free(bitmap.bytes);
    return EXIT_FAILURE;
  }
  SET_THREAD_ARGS(
      args,
      0,
      bitmap,
      batch,
      reals,
      imags,
      row_size,
      img_size,
      top_mask,
      masks
  );
  compute(args);
#else
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
	  batch,
	  reals,
	  imags,
	  row_size,
	  eq_parts + 1,
	  top_mask,
	  masks
      );
      offset += eq_parts + 1;
      if (thrd_create(threads + t, compute, args) != thrd_success) {
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
	batch,
	reals,
	imags,
	row_size,
	eq_parts,
	top_mask,
	masks
    );
    offset += eq_parts;
    if (thrd_create(threads + t, compute, args) != thrd_success) {
      goto thread_fail;
    }
  }
  for (t = 0; t < N_THREADS; ++t) {
    if (thrd_join(threads[t], NULL) != thrd_success) {
      PRINTF_ERR("Thread #%zu failed, exiting.", t + 1);
      failed = true;
    }
  }
#endif
  free(imags);
  free(reals);
#ifndef NO_THREADS
  if (failed) {
    goto late_fail;
  }
#endif
  if (
     printf("P4\n%lu %lu\n", img_size, img_size) < 0 ||
     (
      batch_long ?
      fwrite(bitmap.words, sizeof(unsigned long), bitmap_size, stdout) :
      fwrite(bitmap.bytes, sizeof(unsigned char), bitmap_size, stdout)
     ) != bitmap_size
  ) {
    PUT_ERR("Data output failed, exiting.");
    goto late_fail;
  }
  batch_long ? free(bitmap.words) : free(bitmap.bytes);
  return EXIT_SUCCESS;
#ifndef NO_THREADS
 thread_fail:
  /* Uh oh, let previous threads finish and exit */
  PUT_ERR("Thread creation failed (out of memory?)");
  if (t) {
    --t;
    for (; t < SIZE_MAX; --t) {
      /* Don't care about thread status because nothing will be written */
      thrd_join(threads[t], NULL);
    }
  }
  free(imags);
  free(reals);
#endif
 late_fail:
  batch_long ? free(bitmap.words) : free(bitmap.bytes);
  return EXIT_FAILURE;
}

#undef PUT_ERRNO
#undef PRINTF_ERR
#undef PUT_WARN
#undef PUT_ERR
#undef SET_THREAD_ARGS
#undef MAX_IMG_SIZE
#undef EXPAND
#undef STRINGIFY
#undef UNREACHABLE
#undef NO_THREADS
#undef N_THREADS
