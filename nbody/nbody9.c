/*
 * nbody.c - N-body simulation
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

/* This is a slightly modified translation of C gcc #9. */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 202311L
#error "This program requires C23 support."
#endif
#define OP_4X4(dest, src1, src2, eq, op)	\
  dest 0] eq src1 0] op src2 0];		\
  dest 1] eq src1 1] op src2 1];		\
  dest 2] eq src1 2] op src2 2];		\
  dest 3] eq src1 3] op src2 3]
#define OP_4X1(dest, src1, src2, eq, op)	\
  dest 0] eq src1 0] op src2;		\
  dest 1] eq src1 1] op src2;		\
  dest 2] eq src1 2] op src2;		\
  dest 3] eq src1 3] op src2
#define HADD_4X4(dest, src1, src2)	  \
  dest[0] = src1[0] + src1[1];		  \
  dest[1] = src2[0] + src2[1];		  \
  dest[2] = src1[2] + src1[3];		  \
  dest[3] = src2[2] + src2[3]
#define BLEND_SUM_4X4(dest, src1, src2)	\
  dest 0] = src1[0] + src1[2];		\
  dest 1] = src1[1] + src1[3];		\
  dest 2] = src2[0] + src2[2];		\
  dest 3] = src2[1] + src2[3]
#define COPY_4(dest, src) \
  dest 0] = src 0];	  \
  dest 1] = src 1];	  \
  dest 2] = src 2];	  \
  dest 3] = src 3]

enum body_consts {
  N_BODIES = 5,
  N_COMBNS = (N_BODIES - 1) * N_BODIES / 2,
  VEC_SIZE = 4
};

static inline void rsqrt_vec(
  double vec[static restrict const VEC_SIZE]
) {
  double const tmp[VEC_SIZE] = {
    1.f / sqrtf((float)vec[0]),
    1.f / sqrtf((float)vec[1]),
    1.f / sqrtf((float)vec[2]),
    1.f / sqrtf((float)vec[3]),
  };
  vec[0] = tmp[0] * (1.5 - 0.5 * vec[0] * tmp[0] * tmp[0]);
  vec[1] = tmp[1] * (1.5 - 0.5 * vec[1] * tmp[1] * tmp[1]);
  vec[2] = tmp[2] * (1.5 - 0.5 * vec[2] * tmp[2] * tmp[2]);
  vec[3] = tmp[3] * (1.5 - 0.5 * vec[3] * tmp[3] * tmp[3]);
}

static inline void rsqrt_dist(
  double res[static restrict const N_COMBNS + 2],
  double dst[static restrict const N_COMBNS + 2][VEC_SIZE],
  double pos[static restrict const N_BODIES][VEC_SIZE]
) {
  int i, j, k;
  double tmp1[4][VEC_SIZE], tmp2[VEC_SIZE], tmp3[VEC_SIZE];
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      OP_4X4(dst[k][, pos[i][, pos[j][, =, -);
    }
  }
  for (k = 0; k < N_COMBNS; k += 4) {
    OP_4X4(tmp1[0][, dst[k][, dst[k][, =, *);
    OP_4X4(tmp1[1][, dst[k + 1][, dst[k + 1][, =, *);
    OP_4X4(tmp1[2][, dst[k + 2][, dst[k + 2][, =, *);
    OP_4X4(tmp1[3][, dst[k + 3][, dst[k + 3][, =, *);
    HADD_4X4(tmp2, tmp1[0], tmp1[1]);
    HADD_4X4(tmp3, tmp1[2], tmp1[3]);
    BLEND_SUM_4X4(tmp1[0][, tmp2, tmp3);
    rsqrt_vec(tmp1[0]);
    COPY_4(res[k + , tmp1[0][);
  }
}

static double energy(
  double pos[static restrict const N_BODIES][VEC_SIZE],
  double vel[static restrict const N_BODIES][VEC_SIZE],
  double const mass[static restrict const N_BODIES]
) {
  auto e = .0;
  double dists[N_COMBNS + 2][VEC_SIZE] = {
    [N_BODIES] = {},
    [N_BODIES + 1] = {},
    [N_BODIES + 2] = {},
    [N_COMBNS] = {1., 1., 1., 1.},
    [N_COMBNS + 1] = {1., 1., 1., 1.}
  };
  double mags[N_COMBNS + 2], tmp1[VEC_SIZE], tmp2[VEC_SIZE];
  int i, j, k;
  for (i = 0; i < N_BODIES; ++i) {
    OP_4X4(dists[i][, vel[i][, vel[i][, =, *);
  }
  for (i = 0; i < N_BODIES; i += 4) {
    HADD_4X4(tmp1, dists[i], dists[i + 1]);
    HADD_4X4(tmp2, dists[i + 2], dists[i + 3]);
    BLEND_SUM_4X4(mags[i +, tmp1, tmp2);
  }
  for (i = 0; i < N_BODIES; ++i) {
    e += 0.5 * mass[i] * mags[i];
  }
  rsqrt_dist(mags, dists, pos);
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      e -= mass[i] * mass[j] * mags[k];
    }
  }
  return e;
}

static inline void advance(
  long const n,
  double const dt,
  double const mass[static restrict const N_BODIES],
  double pos[static restrict const N_BODIES][VEC_SIZE],
  double vel[static restrict const N_BODIES][VEC_SIZE]
) {
  double const ticks[VEC_SIZE] = {dt, dt, dt, dt};
  double dists[N_COMBNS + 2][VEC_SIZE] = {
    [N_COMBNS] = {1., 1., 1., 1.},
    [N_COMBNS + 1] = {1., 1., 1., 1.}
  };
  double mags[N_COMBNS + 2], tmp1[VEC_SIZE], tmp2[VEC_SIZE];
  long steps;
  int i, j, k;
  for (steps = 0; steps < n; ++steps) {
    rsqrt_dist(mags, dists, pos);
    for (i = 0; i < N_COMBNS; i += 4) {
      OP_4X4(tmp1[, mags[i +, mags[i +, =, *);
      OP_4X4(tmp2[, mags[i +, ticks[, =, *);
      OP_4X4(mags[i +, tmp1[, tmp2[, =, *);
    }
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	OP_4X1(tmp1[, dists[k][, mags[k], =, *);
        COPY_4(tmp2[, tmp1[);
        OP_4X1(vel[i][, tmp1[, mass[j], -=, *);
        OP_4X1(vel[j][, tmp2[, mass[i], +=, *);
      }
    }
    for (i = 0; i < N_BODIES; ++i) {
      OP_4X4(pos[i][, vel[i][, ticks[, +=, *);
    }
  }
}

int main(int argc, char *argv[]) {
  constexpr auto pi = 3.141'592'653'589'793;
  constexpr auto solar_mass = 4. * pi * pi;
  constexpr auto days_per_year = 365.24;
  constexpr double mass[N_BODIES] = {
    solar_mass,                                   /*     Sun */
    0.000'954'791'938'424'326'609 * solar_mass,   /* Jupiter */
    0.000'285'885'980'666'130'812 * solar_mass,   /*  Saturn */
    0.000'043'662'440'433'515'629'8 * solar_mass, /*  Uranus */
    0.000'051'513'890'204'661'145'1 * solar_mass  /* Neptune */
  };
  char *end = nullptr;
  double pos[N_BODIES][VEC_SIZE] = {
    /* Sun */
    {},
    /* Jupiter */
    {
      .0,
      4.841'431'442'464'720'9,
      -1.160'320'044'027'428'39,
      -0.103'622'044'471'123'109
    },
    /* Saturn */
    {
      .0,
      8.343'366'718'244'579'87,
      4.124'798'564'124'304'79,
      -0.403'523'417'114'321'381
    },
    /* Uranus */
    {
      .0,
      12.894'369'562'139'131,
      -15.111'151'401'698'631'2,
      -0.223'307'578'892'655'734
    },
    /* Neptune */
    {
      .0,
      15.379'697'114'850'916'5,
      -25.919'314'609'987'964'1,
      0.179'258'772'950'371'181
    }
  };
  double vel[N_BODIES][VEC_SIZE] = {
    /* Sun */
    {},
    /* Jupiter */
    {
      .0,
      0.001'660'076'642'744'036'94 * days_per_year,
      0.007'699'011'184'197'404'25 * days_per_year,
      -0.000'069'046'001'697'206'302'3 * days_per_year
    },
    /* Saturn */
    {
      .0,
      -0.002'767'425'107'268'624'11 * days_per_year,
      0.004'998'528'012'349'172'38 * days_per_year,
      0.000'023'041'729'757'376'392'9 * days_per_year
    },
    /* Uranus */
    {
      .0,
      0.002'964'601'375'647'616'18 * days_per_year,
      0.002'378'471'739'594'809'5 * days_per_year,
      -0.000'029'658'956'854'023'755'6 * days_per_year
    },
    /* Neptune */
    {
      .0,
      0.002'680'677'724'903'893'22 * days_per_year,
      0.001'628'241'700'382'422'95 * days_per_year,
      -0.000'095'159'225'451'971'587 * days_per_year
    }
  };
  double e;
  long n;
  {
    double tmp[4] = {};
    for (int i = 1; i < N_BODIES; ++i) {
      OP_4X1(tmp[, vel[i][, mass[i], +=, *);
    }
    OP_4X1(vel[0][, -tmp[, mass[0], =, /);
  }
  e = energy(pos, vel, mass);
  if (argc != 2) {
    puts("Usage: ./nbody iter_count");
    return EXIT_FAILURE;
  } else if ((n = strtol(argv[1], &end, 10)) == 0) {
    puts("\"iter_count\" must be a positive integer");
    return EXIT_FAILURE;
  } else if (errno == ERANGE) {
    perror("Error");
    return EXIT_FAILURE;
  } else {
    advance(n, 0.01, mass, pos, vel);
    printf("%.9f\n%.9f\n", e, energy(pos, vel, mass));
    return EXIT_SUCCESS;
  }
}
