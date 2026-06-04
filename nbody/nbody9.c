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
  double tmp[4][VEC_SIZE];
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      dst[k][0] = pos[i][0] - pos[j][0];
      dst[k][1] = pos[i][1] - pos[j][1];
      dst[k][2] = pos[i][2] - pos[j][2];
      dst[k][3] = pos[i][3] - pos[j][3];
    }
  }
  for (k = 0; k < N_COMBNS; k += 4) {
    tmp[0][0] = dst[k][0] * dst[k][0];
    tmp[0][1] = dst[k][1] * dst[k][1];
    tmp[0][2] = dst[k][2] * dst[k][2];
    tmp[0][3] = dst[k][3] * dst[k][3];

    tmp[1][0] = dst[k + 1][0] * dst[k + 1][0];
    tmp[1][1] = dst[k + 1][1] * dst[k + 1][1];
    tmp[1][2] = dst[k + 1][2] * dst[k + 1][2];
    tmp[1][3] = dst[k + 1][3] * dst[k + 1][3];

    tmp[2][0] = dst[k + 2][0] * dst[k + 2][0];
    tmp[2][1] = dst[k + 2][1] * dst[k + 2][1];
    tmp[2][2] = dst[k + 2][2] * dst[k + 2][2];
    tmp[2][3] = dst[k + 2][3] * dst[k + 2][3];

    tmp[3][0] = dst[k + 3][0] * dst[k + 3][0];
    tmp[3][1] = dst[k + 3][1] * dst[k + 3][1];
    tmp[3][2] = dst[k + 3][2] * dst[k + 3][2];
    tmp[3][3] = dst[k + 3][3] * dst[k + 3][3];
    
    tmp[0][0] = tmp[0][0] + tmp[0][1];
    tmp[0][1] = tmp[1][0] + tmp[1][1];
    tmp[0][2] = tmp[0][2] + tmp[0][3];
    tmp[0][3] = tmp[1][2] + tmp[1][3];

    tmp[2][0] = tmp[2][0] + tmp[2][1];
    tmp[2][1] = tmp[3][0] + tmp[3][1];
    tmp[2][2] = tmp[2][2] + tmp[2][3];
    tmp[2][3] = tmp[3][2] + tmp[3][3];

    tmp[0][0] = tmp[0][0] + tmp[0][2];
    tmp[0][1] = tmp[0][1] + tmp[0][3];
    tmp[0][2] = tmp[2][0] + tmp[2][2];
    tmp[0][3] = tmp[2][1] + tmp[2][3];

    rsqrt_vec(tmp[0]);
    res[k] = tmp[0][0];
    res[k + 1] = tmp[0][1];
    res[k + 2] = tmp[0][2];
    res[k + 3] = tmp[0][3];
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
  double mags[N_COMBNS + 2];
  int i, j, k;
  for (i = 0; i < N_BODIES; ++i) {
    dists[i][0] = vel[i][0] * vel[i][0];
    dists[i][1] = vel[i][1] * vel[i][1];
    dists[i][2] = vel[i][2] * vel[i][2];
    dists[i][3] = vel[i][3] * vel[i][3];
  }
  for (i = 0; i < N_BODIES; i += 4) {
    dists[i][0] = dists[i][0] + dists[i][1];
    dists[i][1] = dists[i + 1][0] + dists[i + 1][1];
    dists[i][2] = dists[i][2] + dists[i][3];
    dists[i][3] = dists[i + 1][2] + dists[i + 1][3];
    dists[i + 2][0] = dists[i + 2][0] + dists[i + 2][1];
    dists[i + 2][1] = dists[i + 3][0] + dists[i + 3][1];
    dists[i + 2][2] = dists[i + 2][2] + dists[i + 2][3];
    dists[i + 2][3] = dists[i + 3][2] + dists[i + 3][3];
    mags[i] = dists[i][0] + dists[i][2];
    mags[i + 1] = dists[i][1] + dists[i][3];
    mags[i + 2] = dists[i + 2][0] + dists[i + 2][2];
    mags[i + 3] = dists[i + 2][1] + dists[i + 2][3];
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
      tmp1[0] = mags[i] * mags[i];
      tmp1[1] = mags[i + 1] * mags[i + 1];
      tmp1[2] = mags[i + 2] * mags[i + 2];
      tmp1[3] = mags[i + 3] * mags[i + 3];
      tmp2[0] = mags[i] * ticks[0];
      tmp2[1] = mags[i + 1] * ticks[1];
      tmp2[2] = mags[i + 2] * ticks[2];
      tmp2[3] = mags[i + 3] * ticks[3];
      mags[i] = tmp1[0] * tmp2[0];
      mags[i + 1] = tmp1[1] * tmp2[1];
      mags[i + 2] = tmp1[2] * tmp2[2];
      mags[i + 3] = tmp1[3] * tmp2[3];
    }
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	tmp1[0] = dists[k][0] * mags[k];
	tmp1[1] = dists[k][1] * mags[k];
	tmp1[2] = dists[k][2] * mags[k];
	tmp1[3] = dists[k][3] * mags[k];
	tmp2[0] = tmp1[0];
	tmp2[1] = tmp1[1];
	tmp2[2] = tmp1[2];
	tmp2[3] = tmp1[3];
	vel[i][0] -= tmp1[0] * mass[j];
	vel[i][1] -= tmp1[1] * mass[j];
	vel[i][2] -= tmp1[2] * mass[j];
	vel[i][3] -= tmp1[3] * mass[j];
	vel[j][0] += tmp2[0] * mass[i];
	vel[j][1] += tmp2[1] * mass[i];
	vel[j][2] += tmp2[2] * mass[i];
	vel[j][3] += tmp2[3] * mass[i];
      }
    }
    for (i = 0; i < N_BODIES; ++i) {
      pos[i][0] += vel[i][0] * ticks[0];
      pos[i][1] += vel[i][1] * ticks[1];
      pos[i][2] += vel[i][2] * ticks[2];
      pos[i][3] += vel[i][3] * ticks[3];
    }
  }
}

int main(int argc, char *argv[]) {
  constexpr auto pi = 3.141592653589793;
  constexpr auto solar_mass = 4. * pi * pi;
  constexpr auto days_per_year = 365.24;
  constexpr double mass[N_BODIES] = {
    solar_mass,                            /*     Sun */
    0.000954791938424326609 * solar_mass,  /* Jupiter */
    0.000285885980666130812 * solar_mass,  /*  Saturn */
    0.0000436624404335156298 * solar_mass, /*  Uranus */
    0.0000515138902046611451 * solar_mass  /* Neptune */
  };
  char *end = nullptr;
  double pos[N_BODIES][VEC_SIZE] = {
    /* Sun */
    {},
    /* Jupiter */
    { .0, 4.8414314424647209, -1.16032004402742839, -0.103622044471123109 },
    /* Saturn */
    { .0, 8.34336671824457987, 4.12479856412430479, -0.403523417114321381 },
    /* Uranus */
    { .0, 12.894369562139131, -15.1111514016986312, -0.223307578892655734 },
    /* Neptune */
    { .0, 15.3796971148509165, -25.9193146099879641, 0.179258772950371181 }
  };
  double vel[N_BODIES][VEC_SIZE] = {
    /* Sun */
    {},
    /* Jupiter */
    {
      .0,
      0.00166007664274403694 * days_per_year,
      0.00769901118419740425 * days_per_year,
      -0.0000690460016972063023 * days_per_year
    },
    /* Saturn */
    {
      .0,
      -0.00276742510726862411 * days_per_year,
      0.00499852801234917238 * days_per_year,
      0.0000230417297573763929 * days_per_year
    },
    /* Uranus */
    {
      .0,
      0.00296460137564761618 * days_per_year,
      0.0023784717395948095 * days_per_year,
      -0.0000296589568540237556 * days_per_year
    },
    /* Neptune */
    {
      .0,
      0.00268067772490389322 * days_per_year,
      0.00162824170038242295 * days_per_year,
      -0.000095159225451971587 * days_per_year
    }
  };
  double e;
  long n;
  {
    double tmp[4] = {};
    for (int i = 1; i < N_BODIES; ++i) {
      tmp[0] += vel[i][0] * mass[i];
      tmp[1] += vel[i][1] * mass[i];
      tmp[2] += vel[i][2] * mass[i];
      tmp[3] += vel[i][3] * mass[i];
    }
    vel[0][0] = -tmp[0] / mass[0];
    vel[0][1] = -tmp[1] / mass[0];
    vel[0][2] = -tmp[2] / mass[0];
    vel[0][3] = -tmp[3] / mass[0];
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
