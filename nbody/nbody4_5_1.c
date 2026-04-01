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

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum body_consts {
  N_BODIES = 5,
  N_PARAMS = 7,
  N_COMBS = (N_BODIES - 1) * N_BODIES / 2
};

static inline void offset_momentum(
  double (* restrict const bodies)[N_PARAMS]
) {
  int const sun_idx = N_BODIES - 1;
  int i;
  double mass, px = 0., py = 0., pz = 0.;
  for (i = 0; i < sun_idx; ++i) {
    mass = bodies[i][6];
    px += bodies[i][3] * mass;
    py += bodies[i][4] * mass;
    pz += bodies[i][5] * mass;
  }
  mass = bodies[sun_idx][6];
  bodies[sun_idx][3] = -px / mass;
  bodies[sun_idx][4] = -py / mass;
  bodies[sun_idx][5] = -pz / mass;
}

/*
 * Fast inverse square root
 * Why: sqrt (the C math function) is impure, so it cannot be optimized
 * out, and this prevents the first call to "energy" from being
 * calculated at compile time.
 */
static inline double pure_rsqrt(double x) {
  union { double f; unsigned long long i; } bits = { .f = x };
  bits.i = 0x5fe6eb50c7b537a9ull - (bits.i >> 1);
  bits.f *= 1.5 - 0.5 * x * bits.f * bits.f;
  bits.f *= 1.5 - 0.5 * x * bits.f * bits.f;
  return bits.f * (1.5 - 0.5 * x * bits.f * bits.f);
}

static inline double energy(
  double (* restrict const bodies)[N_PARAMS],
  double * restrict const mags
) {
  int i, j, k;
  double e = 0., dx, dy, dz, d2, dist, vx, vy, vz;
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      dx = bodies[i][0] - bodies[j][0];
      dy = bodies[i][1] - bodies[j][1];
      dz = bodies[i][2] - bodies[j][2];
      mags[k] = dx * dx + dy * dy + dz * dz;
    }
  }
  for (k = 0; k < N_COMBS; ++k) {
    d2 = mags[k];
    dist = 1.f / sqrtf((float)d2);
    mags[k] = dist * (1.5 - 0.5 * d2 * dist * dist);
  }
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      e -= bodies[i][6] * bodies[j][6] * mags[k];
    }
    vx = bodies[i][3];
    vy = bodies[i][4];
    vz = bodies[i][5];
    e += 0.5 * bodies[i][6] * (vx * vx + vy * vy + vz * vz);
  }
  return e;
}

static inline void advance(
  double (* restrict const bodies)[N_PARAMS],
  double (* restrict const dists)[3],
  double * restrict const mags,
  double const dt
) {
  int i, j, k;
  double dx, dy, dz, d2, dist, mag, mm;
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      dists[k][0] = bodies[i][0] - bodies[j][0];
      dists[k][1] = bodies[i][1] - bodies[j][1];
      dists[k][2] = bodies[i][2] - bodies[j][2];
    }
  }
  for (k = 0; k < N_COMBS; ++k) {
    dx = dists[k][0];
    dy = dists[k][1];
    dz = dists[k][2];
    d2 = dx * dx + dy * dy + dz * dz;
    dist = 1.f / sqrtf((float)d2);
    dist *= 1.5 - 0.5 * d2 * dist * dist;
    mags[k] = dt * (dist * dist * dist);
  }
  for (i = 0, k = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j, ++k) {
      mag = mags[k];
      dx = dists[k][0];
      dy = dists[k][1];
      dz = dists[k][2];
      mm = bodies[j][6] * mag;
      bodies[i][3] -= dx * mm;
      bodies[i][4] -= dy * mm;
      bodies[i][5] -= dz * mm;
      mm = bodies[i][6] * mag;
      bodies[j][3] += dx * mm;
      bodies[j][4] += dy * mm;
      bodies[j][5] += dz * mm;
    }
    bodies[i][0] += bodies[i][3] * dt;
    bodies[i][1] += bodies[i][4] * dt;
    bodies[i][2] += bodies[i][5] * dt;
  }
}

int main(int argc, char *argv[]) {
  double const pi = 3.141592653589793;
  double const solar_mass = 4. * pi * pi;
  double const days_per_year = 365.24;
  char *end = NULL;
  double
    (* const bodies)[N_PARAMS] = malloc((size_t)N_BODIES * sizeof *bodies),
    (* const dists)[3] = malloc((size_t)N_COMBS * sizeof *dists),
    e1 = 0.;
  unsigned long iterations, step;
  double e2, mags[N_COMBS];
  if (!bodies || !dists) {
    puts("FAIL: could not allocate matrices, not enough memory.");
    goto fail;
  }
  if (argc < 2) {
    puts("FAIL: Please provide the amount of iterations.");
    goto fail;
  }
  errno = 0;
  iterations = strtoul(argv[1], &end, 10);
  if (errno == ERANGE) {
    perror("FAIL");
    goto fail;
  }
  if (!iterations) {
    puts("FAIL: Zero iterations (is the argument malformed?)");
    goto fail;
  }
  /*
   * Astronomical bodies:
   *  - doubles 0, 1 and 2 form the position
   *  - doubles 3, 4 and 5 form the velocity
   *  - double 6 represents the mass
   */
  /* Jupiter */
  bodies[0][0] = 4.8414314424647209;
  bodies[0][1] = -1.16032004402742839;
  bodies[0][2] = -0.103622044471123109;
  bodies[0][3] = 0.00166007664274403694 * days_per_year;
  bodies[0][4] = 0.00769901118419740425 * days_per_year;
  bodies[0][5] = -0.0000690460016972063023 * days_per_year;
  bodies[0][6] = 0.000954791938424326609 * solar_mass;
  /* Saturn */
  bodies[1][0] = 8.34336671824457987;
  bodies[1][1] = 4.12479856412430479;
  bodies[1][2] = -0.403523417114321381;
  bodies[1][3] = -0.00276742510726862411 * days_per_year;
  bodies[1][4] = 0.00499852801234917238 * days_per_year;
  bodies[1][5] = 0.0000230417297573763929 * days_per_year;
  bodies[1][6] = 0.000285885980666130812 * solar_mass;
  /* Uranus */
  bodies[2][0] = 12.894369562139131;
  bodies[2][1] = -15.1111514016986312;
  bodies[2][2] = -0.223307578892655734;
  bodies[2][3] = 0.00296460137564761618 * days_per_year;
  bodies[2][4] = 0.0023784717395948095 * days_per_year;
  bodies[2][5] = -0.0000296589568540237556 * days_per_year;
  bodies[2][6] = 0.0000436624404335156298 * solar_mass;
  /* Neptune */
  bodies[3][0] = 15.3796971148509165;
  bodies[3][1] = -25.9193146099879641;
  bodies[3][2] = 0.179258772950371181;
  bodies[3][3] = 0.00268067772490389322 * days_per_year;
  bodies[3][4] = 0.00162824170038242295 * days_per_year;
  bodies[3][5] = -0.000095159225451971587 * days_per_year;
  bodies[3][6] = 0.0000515138902046611451 * solar_mass;
  /* Sun */
  bodies[4][0] = 0.;
  bodies[4][1] = 0.;
  bodies[4][2] = 0.;
  bodies[4][3] = 0.;
  bodies[4][4] = 0.;
  bodies[4][5] = 0.;
  bodies[4][6] = solar_mass;
  offset_momentum(bodies);
  /*
   * Manually inline first call to energy and replace sqrt()
   * with the custom reciprocal square root to convince gcc
   * to compute the result at compile time.
   */
  {
    int i, j, k;
    double dx, dy, dz, vx, vy, vz;
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	dx = bodies[i][0] - bodies[j][0];
	dy = bodies[i][1] - bodies[j][1];
	dz = bodies[i][2] - bodies[j][2];
	mags[k] = dx * dx + dy * dy + dz * dz;
      }
    }
    for (k = 0; k < N_COMBS; ++k) {
      mags[k] = pure_rsqrt(mags[k]);
    }
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	e1 -= bodies[i][6] * bodies[j][6] * mags[k];
      }
      vx = bodies[i][3];
      vy = bodies[i][4];
      vz = bodies[i][5];
      e1 += 0.5 * bodies[i][6] * (vx * vx + vy * vy + vz * vz);
    }
  }
  for (step = 0; step < iterations; ++step) {
    advance(bodies, dists, mags, 0.01);
  }
  e2 = energy(bodies, mags);
  printf("%.9f\n%.9f\n", e1, e2);
  free(bodies);
  free(dists);
  return EXIT_SUCCESS;
 fail:
  free(bodies);
  free(dists);
  return EXIT_FAILURE;
}
