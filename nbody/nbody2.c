#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef __GNUC__
#pragma STDC FENV_ACCESS OFF
#pragma STDC FP_CONTRACT ON
#endif

enum body_consts {
  N_BODIES = 5,
  N_PARAMS = 7,
  N_COMBS = (N_BODIES - 1) * N_BODIES / 2
};

static inline void offset_momentum(
  double (* restrict const bodies)[N_BODIES][N_PARAMS]
) {
  int const sun_idx = N_BODIES - 1;
  int i;
  double mass, px = 0., py = 0., pz = 0.;
  for (i = 0; i < sun_idx; ++i) {
    mass = (*bodies)[i][6];
    px += (*bodies)[i][3] * mass;
    py += (*bodies)[i][4] * mass;
    pz += (*bodies)[i][5] * mass;
  }
  mass = (*bodies)[sun_idx][6];
  (*bodies)[sun_idx][3] = -px / mass;
  (*bodies)[sun_idx][4] = -py / mass;
  (*bodies)[sun_idx][5] = -pz / mass;
}

static inline double energy(
  double (* restrict const bodies)[N_BODIES][N_PARAMS]
) {
  int i, j;
  double e = 0., dx, dy, dz, vx, vy, vz;
  for (i = 0; i < N_BODIES; ++i) {
    for (j = i + 1; j < N_BODIES; ++j) {
      dx = (*bodies)[i][0] - (*bodies)[j][0];
      dy = (*bodies)[i][1] - (*bodies)[j][1];
      dz = (*bodies)[i][2] - (*bodies)[j][2];
      e -= (*bodies)[i][6] * (*bodies)[j][6] / sqrt(dx * dx + dy * dy + dz * dz);
    }
    vx = (*bodies)[i][3];
    vy = (*bodies)[i][4];
    vz = (*bodies)[i][5];
    e += 0.5 * (*bodies)[i][6] * (vx * vx + vy * vy + vz * vz);
  }
  return e;
}

static inline void advance(
  double (* restrict const bodies)[N_BODIES][N_PARAMS],
  double const dt,
  unsigned long const steps
) {
  double
    (* const dists)[3] = malloc((size_t)N_COMBS * sizeof *dists),
    * const mags = malloc((size_t)N_COMBS * sizeof *mags);
  int i, j, k;
  double d2, dist, mm;
  unsigned long step;
  if (!dists || !mags) {
    /* If the program is out of memory, logging most likely won't work */
    abort();
  }
  for (step = 0ul; step < steps; ++step) {
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	dists[k][0] = (*bodies)[i][0] - (*bodies)[j][0];
	dists[k][1] = (*bodies)[i][1] - (*bodies)[j][1];
	dists[k][2] = (*bodies)[i][2] - (*bodies)[j][2];
      }
    }
    for (k = 0; k < N_COMBS; ++k) {
      d2 =
	dists[k][0] * dists[k][0] +
	dists[k][1] * dists[k][1] +
	dists[k][2] * dists[k][2];
      dist = 1.f / sqrtf((float)d2);
      dist *= 1.5 - 0.5 * d2 * dist * dist;
      mags[k] = dt * (dist * dist * dist);
    }
    for (i = 0, k = 0; i < N_BODIES; ++i) {
      for (j = i + 1; j < N_BODIES; ++j, ++k) {
	mm = (*bodies)[j][6] * mags[k];
	(*bodies)[i][3] -= dists[k][0] * mm;
	(*bodies)[i][4] -= dists[k][1] * mm;
	(*bodies)[i][5] -= dists[k][2] * mm;
	mm = (*bodies)[i][6] * mags[k];
	(*bodies)[j][3] += dists[k][0] * mm;
	(*bodies)[j][4] += dists[k][1] * mm;
	(*bodies)[j][5] += dists[k][2] * mm;
      }
      (*bodies)[i][0] += (*bodies)[i][3] * dt;
      (*bodies)[i][1] += (*bodies)[i][4] * dt;
      (*bodies)[i][2] += (*bodies)[i][5] * dt;
    }
  }
  free(dists);
  free(mags);
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
  double bodies[N_BODIES][N_PARAMS] = {
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
  };
  char *end = NULL;
  unsigned long iterations;
  double e1, e2;
  if (argc < 2) {
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
  offset_momentum(&bodies);
  e1 = energy(&bodies);
  advance(&bodies, 0.01, iterations);
  e2 = energy(&bodies);
  printf("%.9f\n%.9f\n", e1, e2);
  return EXIT_SUCCESS;
}
