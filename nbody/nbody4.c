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
  double bodies[static restrict const N_BODIES][N_PARAMS]
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

static inline double energy(
  double bodies[static restrict const N_BODIES][N_PARAMS],
  double mags[static restrict const N_COMBS]
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
  double bodies[static restrict const N_BODIES][N_PARAMS],
  double dists[static restrict const N_COMBS][3],
  double mags[static restrict const N_COMBS],
  double const dt,
  unsigned long const steps
) {
  int i, j, k;
  double dx, dy, dz, d2, dist, mag, mm;
  unsigned long step;
  for (step = 0ul; step < steps; ++step) {
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
}

int main(int argc, char *argv[]) {
  double const pi = 3.141592653589793;
  double const solar_mass = 4. * pi * pi;
  double const days_per_year = 365.24;
  /*
   * Astronomical bodies:
   *  - doubles 0, 1 and 2 form the position
   *  - doubles 3, 4 and 5 form the velocity
   *  - double 6 represents the mass
   */
  double
    bodies[N_BODIES][N_PARAMS] = {
    /* Jupiter */
    {
      4.8414314424647209,
      -1.16032004402742839,
      -0.103622044471123109,
      0.00166007664274403694 * days_per_year,
      0.00769901118419740425 * days_per_year,
      -0.0000690460016972063023 * days_per_year,
      0.000954791938424326609 * solar_mass
    },
    /* Saturn */
    {
      8.34336671824457987,
      4.12479856412430479,
      -0.403523417114321381,
      -0.00276742510726862411 * days_per_year,
      0.00499852801234917238 * days_per_year,
      0.0000230417297573763929 * days_per_year,
      0.000285885980666130812 * solar_mass
    },
    /* Uranus */
    {
      12.894369562139131,
      -15.1111514016986312,
      -0.223307578892655734,
      0.00296460137564761618 * days_per_year,
      0.0023784717395948095 * days_per_year,
      -0.0000296589568540237556 * days_per_year,
      0.0000436624404335156298 * solar_mass
    },
    /* Neptune */
    {
      15.3796971148509165,
      -25.9193146099879641,
      0.179258772950371181,
      0.00268067772490389322 * days_per_year,
      0.00162824170038242295 * days_per_year,
      -0.000095159225451971587 * days_per_year,
      0.0000515138902046611451 * solar_mass
    },
    /* Sun */
    {0., 0., 0., 0., 0., 0., solar_mass}
  },
    (* const dists)[3] = malloc((size_t)N_COMBS * sizeof *dists),
    * const mags = malloc((size_t)N_COMBS * sizeof *mags);
  char *end = NULL;
  unsigned long iterations;
  double e1, e2;
  if (!dists || !mags) {
    puts("FAIL: Cannot allocate memory for arrays.");
    return EXIT_FAILURE;
  } else if (argc < 2) {
    puts("FAIL: Please provide the amount of iterations.");
    return EXIT_FAILURE;
  }
  errno = 0;
  iterations = strtoul(argv[1], &end, 10);
  if (errno == ERANGE) {
    perror("FAIL");
    return EXIT_FAILURE;
  } else if (!iterations) {
    puts("FAIL: Zero iterations (is the argument malformed?)");
    return EXIT_FAILURE;
  }
  offset_momentum(bodies);
  e1 = energy(bodies, mags);
  advance(bodies, dists, mags, 0.01, iterations);
  e2 = energy(bodies, mags);
  printf("%.9f\n%.9f\n", e1, e2);
  free(dists);
  free(mags);
  return EXIT_SUCCESS;
}
