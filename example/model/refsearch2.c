/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <sc_options.h>
#include "model.h"

static int max_ref_level = 3;

typedef struct triangle
{
  double v0[3];
  double v1[3];
  double v2[3];
}
triangle_t;

#if 0
static int
triangulation_is_vertex_inside_aabb (const double * aabb, const double * v)
{
  return ((v[0] >= aabb[0] && v[0] <= aabb[3]) &&
          (v[1] >= aabb[1] && v[1] <= aabb[4]) &&
          (v[2] >= aabb[2] && v[2] <= aabb[5]));
}

static int
triangulation_intersect_model (p4est_topidx_t which_tree,
                               const double aabb[6], void *primitive)
{
  triangle_t * t = (triangle_t *) primitive;

  int is_v0_in = triangulation_is_vertex_inside_aabb (aabb, t->v0);
  int is_v1_in = triangulation_is_vertex_inside_aabb (aabb, t->v1);
  int is_v2_in = triangulation_is_vertex_inside_aabb (aabb, t->v2);

  if ((is_v0_in && is_v1_in && is_v2_in)
   || (!is_v0_in && !is_v1_in && !is_v2_in))
   {
     return 0;
   }
   else {
     return 1;
   }
}

static void
triangulation_setup_model (p4est_model_t ** m)
{
  p4est_model_t  *model = P4EST_ALLOC_ZERO (p4est_model_t, 1);
  model->output_prefix = "triangulation";
  model->conn = p4est_connectivity_new_unitsquare ();
  model->geom = NULL;
  model->intersect = triangulation_intersect_model;
}
#endif

static int
usagerrf (sc_options_t * opt, const char *fmt, ...)
{
  va_list             ap;
  char                msg[BUFSIZ];

  va_start (ap, fmt);
  vsnprintf (msg, BUFSIZ, fmt, ap);
  va_end (ap);

  P4EST_GLOBAL_LERROR ("ERROR/\n");
  P4EST_GLOBAL_LERRORF ("ERROR: %s\n", msg);
  P4EST_GLOBAL_LERROR ("ERROR\\\n");
  return 1;
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  return usagerrf (opt, "%s", msg);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 ue, fa;
  sc_options_t       *opt;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize global context */
  mpicomm = sc_MPI_COMM_WORLD;

  /* set global logging options for p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global application state */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'L', "maxlevel", &max_ref_level, P4EST_QMAXLEVEL,
                      "Maximum refinement level");

    /* proceed in run-once loop for cleaner error checking */
  ue = 0;
  do {
    /* parse command line and assign configuration variables */
    fa = sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);
    if (fa < 0 || fa != argc) {
      ue = usagerr (opt, "invalid option format or non-option argument");
      break;
    }
    P4EST_GLOBAL_PRODUCTIONF ("Manifold dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

    /* check consistency of parameters */
    if (max_ref_level < 1 || max_ref_level > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "maxlevel not between 1 and %d",
                     P4EST_QMAXLEVEL);
    }
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  return 0;
}