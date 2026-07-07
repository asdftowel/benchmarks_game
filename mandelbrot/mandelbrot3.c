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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#if !defined(NAN) || !defined(INFINITY)
#error "This program requires support for NaN and Infinity."
#endif

#ifndef N_THREADS
#define N_THREADS 4 /* Adjust to number of hardware threads */
#elif N_THREADS < 1
#error "N_THREADS must be at least 1."
#endif

#if !defined(__STDC_VERSION__) || \
  __STDC_VERSION__ < 201112L   || \
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
#include <stddef.h>
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

static void calc_sum(
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
  bool result = true;
  int i;
  for (i = 0; i < CHAR_BIT; ++i) {
    if (vec[i] <= 4.) {
      result = false;
      break;
    }
  }
  return result;
}

enum bit_consts {
  ULONG_BITS = CHAR_BIT * sizeof(unsigned long),
  HIGH_SHIFT = ULONG_BITS - CHAR_BIT,
  TOP_BIT = CHAR_BIT - 1
};

static unsigned int mand_char(
    double const * restrict const init_r,
    double const init_i,
    unsigned int const * restrict const masks
) {
  double reals[CHAR_BIT], imags[CHAR_BIT], sums[CHAR_BIT];
  unsigned int result = UCHAR_MAX;
  int i, j;

  for (i = 0; i < CHAR_BIT; ++i) {
    reals[i] = init_r[i];
    imags[i] = init_i;
  }
  
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
  j = TOP_BIT;
  for (i = 0; i < CHAR_BIT; ++i) {
    result &= masks[i] | ((unsigned int)(sums[i] <= 4.) << j);
    --j;
  }
 done:
  return result;
}

static inline unsigned long mand_long(
    double const * restrict const init_r,
    double const init_i,
    unsigned int const * restrict const masks
) {
  unsigned long result = 0;
  int offset = 0;
  size_t i;
  for (i = 0; i < sizeof(unsigned long); ++i) {
    result =
      (result >> CHAR_BIT) |
      (
       (unsigned long)mand_char(init_r + offset, init_i, masks)
       << HIGH_SHIFT
      );
    offset += CHAR_BIT;
  }
  return result;
}

enum batch_size { BYTE, WORD };

union bitmap_type { unsigned char * bytes; unsigned long * words; };

struct thread_args {
  union bitmap_type bitmap;
  enum batch_size batch;
  double const * init_r;
  double const * init_i;
  unsigned long width;
  unsigned long nrows;
  unsigned int const * masks;
};

enum mand_state {
  SUCCEEDED = 0,
  ALLOC_FAILED,
  JOIN_FAILED
};

#ifdef NO_THREADS
static inline int mand_compute(
    union bitmap_type const bitmap,
    enum batch_size const batch,
    unsigned long const nrows,
    unsigned long const width,
    unsigned int const * restrict const masks,
    double const * restrict const init_r,
    double const * restrict const init_i
) {
#else
static int mand_compute(void *args) {
  struct thread_args const * const locals = args;
  union bitmap_type const bitmap = locals->bitmap;
  enum batch_size const batch = locals->batch;
  double const * const init_r = locals->init_r,
    * const init_i = locals->init_i;
  unsigned long const width = locals->width,
    nrows = locals->nrows;
  unsigned int const * const masks = locals->masks;
#endif
  unsigned long row, col, ir_offset;
  size_t row_offset = 0;
  switch (batch) {
  case BYTE:
    for (row = 0; row < nrows; ++row) {
      ir_offset = 0;
      for (col = 0; col < width; ++col) {
        bitmap.bytes[row_offset + col] =
          (unsigned char)mand_char(init_r + ir_offset, init_i[row], masks);
        ir_offset += CHAR_BIT;
      }
      row_offset += width;
    }
    break;
  case WORD:
    for (row = 0; row < nrows; ++row) {
      ir_offset = 0;
      for (col = 0; col < width; ++col) {
        bitmap.words[row_offset + col] =
          mand_long(init_r + ir_offset, init_i[row], masks);
        ir_offset += ULONG_BITS;
      }
      row_offset += width;
    }
    break;
  default:
    UNREACHABLE;
  }
#ifndef NO_THREADS
  free(args);
#endif
  return 0;
}

#ifndef NO_THREADS
#define SET_THREAD_ARGS(name, off, bm, s, r, i, w, rows, m)     \
  switch (s) {                                                  \
  case BYTE:                                                    \
    name->bitmap.bytes = bm.bytes + off * w;                    \
    break;                                                      \
  case WORD:                                                    \
    name->bitmap.words = bm.words + off * w;                    \
    break;                                                      \
  default:                                                      \
    UNREACHABLE;                                                \
  }                                                             \
  name->batch = s;                                              \
  name->init_r = r;                                             \
  name->init_i = i + off;                                       \
  name->width = w;                                              \
  name->nrows = rows;                                           \
  name->masks = m

static inline enum mand_state mand_manage(
    union bitmap_type const bitmap,
    enum batch_size const batch,
    unsigned long const img_size,
    unsigned long const row_size,
    unsigned int const * restrict const masks,
    double const * restrict const reals,
    double const * restrict const imags
) {
  unsigned long eql_parts, remaining;
  size_t img_offset = 0, t;
  thrd_t threads[N_THREADS];
  struct thread_args * args;
  bool status = SUCCEEDED;
  eql_parts = img_size / N_THREADS;
  remaining = img_size % N_THREADS;
  if (!eql_parts) {
    /*
     * If there's not enough work for all threads,
     * it can probably be done by two.
     */
    eql_parts = img_size / 2;
    remaining = img_size % 2;
  }
  for (t = 0; t < N_THREADS; ++t) {
    if (!(args = malloc(sizeof * args))) {
      goto alloc_fail;
    }
    SET_THREAD_ARGS(
        args,
        img_offset,
        bitmap,
        batch,
        reals,
        imags,
        row_size,
        remaining ? eql_parts + 1 : eql_parts,
        masks
    );
    img_offset += remaining ? eql_parts + 1 : eql_parts;
    if (thrd_create(threads + t, mand_compute, args) != thrd_success) {
      goto alloc_fail;
    }
    remaining -= remaining != 0;
  }
  for (t = 0; t < N_THREADS; ++t) {
    if (thrd_join(threads[t], NULL) != thrd_success) {
      status = JOIN_FAILED;
    }
  }
  return status;
 alloc_fail:
  status = ALLOC_FAILED;
  while (t) {
    --t;
    thrd_join(threads[t], NULL);
  }
  return status;
}

#undef SET_THREAD_ARGS
#endif

static inline bool mand_attempt_write(
    union bitmap_type const bitmap,
    enum batch_size const batch,
    unsigned long const img_size,
    size_t const bitmap_size
) {
  bool result = false;
  if (printf("P4\n%lu %lu\n", img_size, img_size) < 0) {
    goto fail;
  }
  switch (batch) {
  case BYTE:
    result =
      fwrite(bitmap.bytes, sizeof(unsigned char), bitmap_size, stdout)
      == bitmap_size;
    break;
  case WORD:
    result =
      fwrite(bitmap.words, sizeof(unsigned long), bitmap_size, stdout)
      == bitmap_size;
    break;
  default:
    UNREACHABLE;
  }
 fail:
  return result;
}
 
#ifndef MAX_IMG_SIZE
#define MAX_IMG_SIZE 16000 /* In pixels */
#endif

#define PASTE(tok, suffix) tok ## suffix
#define PASTE_EX(tok, suffix) PASTE(tok, suffix)
#define PUT_ERR(msg)                                            \
  fputs(__FILE__ ", line " EXPAND(__LINE__) ": Error: " msg "\n", stderr)
#define PUT_WARN(msg)                                           \
  fputs(__FILE__ ", line " EXPAND(__LINE__) ": Warning: " msg "\n", stderr)
#define PRINTF_ERR(fmt, ...)                                            \
  fprintf(                                                              \
      stderr,                                                           \
      __FILE__ ", line " EXPAND(__LINE__) ": Error: " fmt "\n",         \
      __VA_ARGS__                                                       \
  )
#define PUT_ERRNO perror(__FILE__ ", line " EXPAND(__LINE__) ": Error")

int main(int argc, char *argv[]) {
  unsigned long img_size, img_rem, row_size, i;
  unsigned int masks[CHAR_BIT];
  int idx, shift, exit_code = EXIT_FAILURE;
  size_t bitmap_size;
  union bitmap_type bitmap;
  char * end = NULL;
  double inv, * reals, * imags;
  bool batch_long;
  enum batch_size batch;

  shift = TOP_BIT;
  for (idx = 0; idx < CHAR_BIT; ++idx) {
    masks[idx] = ~(1u << shift);
    --shift;
  }

#ifdef _WIN32
  /* Windows silently converts \n into \r\n */
  if (_setmode(_fileno(stdout), _O_BINARY) == -1) {
    PUT_ERRNO;
    goto fail_0;
  }
#endif
  
  if (argc != 2) {
    PUT_ERR("Usage: ./mandelbrot <image size>\n");
    goto fail_0;
  } else if (!(img_size = strtoul(argv[1], &end, 10))) {
    PRINTF_ERR("Zero pixels. Is \"%s\" a number?", argv[1]);
    goto fail_0;
  } else if (errno == ERANGE) {
    PUT_ERRNO;
    goto fail_0;
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
  if (!isnormal(inv = 2. / (double)img_size)) {
    PUT_ERR("Cannot represent inv, result will be incorrect.");
    goto fail_0;
  }
  row_size = img_size / CHAR_BIT;
  if (img_size > PASTE_EX(MAX_IMG_SIZE, ul)) {
    PUT_ERR("Image size larger than compile-time limit.");
    goto fail_0;
  }
  
  if (SIZE_MAX / row_size < img_size) {
    PUT_ERR("Input dimensions would cause overflow, try decreasing the size.");
    goto fail_0;
  }
  bitmap_size = (size_t)img_size * row_size;

  if (SIZE_MAX / sizeof(double) < img_size) {
    PUT_ERR("Allocation sizes would cause overflow, try decreasing the size.");
    goto fail_0;
  } else if (!(reals = malloc(img_size * sizeof(double)))) {
    PUT_ERR("Cannot allocate real numbers (out of memory?)");
    goto fail_0;
  } else if (!(imags = malloc(img_size * sizeof(double)))) {
    PUT_ERR("Cannot allocate imaginary numbers (out of memory?)");
    goto fail_1;
  }

  for (i = 0; i < img_size; ++i) {
    reals[i] = inv * (double)i - 1.5;
    imags[i] = inv * (double)i - 1.;
  }

  if (batch_long) {
    batch = WORD;
    row_size /= sizeof(unsigned long);
    if (!(bitmap.words = malloc(bitmap_size))) {
      PUT_ERR("Cannot allocate bitmap (out of memory?)");
      goto fail_2;
    }
    bitmap_size /= sizeof(unsigned long);
  } else {
    batch = BYTE;
    if (!(bitmap.bytes = malloc(bitmap_size))) {
      PUT_ERR("Cannot allocate bitmap (out of memory?)");
      goto fail_2;
    }
  }

  switch (
#ifdef NO_THREADS
      mand_compute(bitmap, batch, img_size, row_size, masks, reals, imags)
#else
      mand_manage(bitmap, batch, img_size, row_size, masks, reals, imags)
#endif
  ) {
  case SUCCEEDED:
    if (mand_attempt_write(bitmap, batch, img_size, bitmap_size)) {
      exit_code = EXIT_SUCCESS;
    } else {
      PUT_ERR("Data output failed, exiting.");
    }
    break;
  case ALLOC_FAILED:
    PUT_ERR("Cannot allocate thread (out of memory?)");
    break;
  case JOIN_FAILED:
    PUT_ERR("At least one thread terminated incorrectly.");
    break;
  default:
    UNREACHABLE;
  }

  batch_long ? free(bitmap.words) : free(bitmap.bytes);
 fail_2:
  free(imags);
 fail_1:
  free(reals);
 fail_0:
  return exit_code;
}

#undef PUT_ERRNO
#undef PRINTF_ERR
#undef PUT_WARN
#undef PUT_ERR
#undef PASTE_EX
#undef PASTE
#undef MAX_IMG_SIZE
#undef EXPAND
#undef STRINGIFY
#undef UNREACHABLE
#undef NO_THREADS
#undef N_THREADS
