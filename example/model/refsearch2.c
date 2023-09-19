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

/* #if 0 */
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

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 * 
 * Copied from p4est_connectivity.c:
 * p4est_connectivity_getline_upper (FILE * stream).
 */
static char        *
triangulation_getline_upper (FILE * stream)
{
  char               *line = P4EST_ALLOC (char, 1024), *linep = line;
  size_t              lenmax = 1024, len = lenmax;
  int                 c;

  if (line == NULL)
    return NULL;

  for (;;) {
    c = fgetc (stream);
    if (c == EOF && linep == line) {
      P4EST_FREE (linep);
      return NULL;
    }
    c = toupper (c);

    if (--len == 0) {
      char               *linen;

      len = lenmax;
      lenmax *= 2;

      linen = P4EST_REALLOC (linep, char, lenmax);
      if (linen == NULL) {
        P4EST_FREE (linep);
        return NULL;
      }

      line = linen + (line - linep);
      linep = linen;
    }
    if ((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

/** This function goes through the file fin that has .off (Object File Format).
 * Object File Format (.off) files are used to represent the geometry of a
 * model by specifying the polygons of the model's surface. The polygons can
 * have any number of vertices.
 * The .off files in the conform to the following standard. OFF files are all
 * ASCII files beginning with the keyword OFF. The next line states the number
 * of vertices, the number of faces, and the number of edges. The number of
 * edges can be safely ignored.
 *
 * The vertices are listed with x, y, z coordinates, written one per line.
 * After the list of vertices, the faces are listed, with one face per line.
 *
 * ****************************************************************************
 * !!!Our model considers triangulated meshes only, so its all faces should be
 * triangles!!!
 * ****************************************************************************
 *
 * For each face, the number of vertices is specified, followed by indices into
 * the list of vertices. See the examples below.
 * Note that earlier versions of the model files had faces with -1 indices into
 * the vertex list. That was due to an error in the conversion program and
 * should be corrected now.
 * OFF
 * numVertices numFaces numEdges
 * x y z
 * x y z
 *
 * ... numVertices like above
 * NVertices v1 v2 v3 ... vN
 * MVertices v1 v2 v3 ... vM
 * ... numFaces like above
 *
 * Note that vertices are numbered starting at 0 (not starting at 1).
 * A simple example for a cube (that is not supported by this model):
 *
 * OFF
 * 8 6 0
 * -0.500000 -0.500000 0.500000
 * 0.500000 -0.500000 0.500000
 * -0.500000 0.500000 0.500000
 * 0.500000 0.500000 0.500000
 * -0.500000 0.500000 -0.500000
 * 0.500000 0.500000 -0.500000
 * -0.500000 -0.500000 -0.500000
 * 0.500000 -0.500000 -0.500000
 * 4 0 1 3 2
 * 4 2 3 5 4
 * 4 4 5 7 6
 * 4 6 7 1 0
 * 4 1 7 5 3
 * 4 6 0 2 4
*/
static int
triangulation_read_off_file_stream (p4est_model_t * m, FILE * fin)
{
  char * line;
  int                 lines_read = 0, v = 0, f = 0;
  int retval;
  size_t num_vertices, num_edges, v2f, vid[3];
  triangle_t * faces;
  double * vertices;
  double x, y, z;
  double min_x, min_y, min_z;
  double max_x, max_y, max_z;
  double axis_scale = 1.;

  for (;;) {
    line = triangulation_getline_upper (fin);

    if (line == NULL) {
      break;
    }

    ++lines_read;

    /* check for file format */
    if (lines_read == 1) {
      if (!strstr (line, "OFF")) {
        P4EST_LERROR ("Wrong file format to read");
        P4EST_FREE (line);
        return 0;
      }
      continue;
    }

    /* check for number of vertices and faces in the object to read */
    if (lines_read == 2) {
      retval
        = sscanf (line, "%lu, %lu, %lu", &num_vertices, &m->num_prim, &num_edges);

      if (retval != 3) {
        P4EST_LERROR ("Wrong file format to read");
        P4EST_FREE (line);
        return 0;
      }
      if (num_vertices < 3) {
        P4EST_LERROR ("Not enough vertices to build a triangle");
        P4EST_FREE (line);
        return 0;
      }
      if (m->num_prim == 0) {
        P4EST_LERROR ("No primitives in the file");
        P4EST_FREE (line);
        return 0;
      }

      /* allocate memory for vertices and face (triangles (primitives)) */
      vertices = P4EST_ALLOC (double, 3 * num_vertices);
      faces = m->primitives = P4EST_ALLOC (triangle_t, m->num_prim);

      continue;
    }

    if (v < num_vertices) {
      /* read vertices */
      retval = sscanf (line, "%lf, %lf, %lf", &x, &y, &z);
      if (retval != 3) {
        P4EST_LERROR ("Premature end of file");
        P4EST_FREE (line);
        return 0;
      }
      if (v == 0) {
        min_x = max_x = x;
        min_y = max_y = y;
        min_z = max_z = z;
      }
      min_x = SC_MIN (min_x, x);
      min_y = SC_MIN (min_y, y);
      min_z = SC_MIN (min_z, z);

      min_x = SC_MAX (max_x, x);
      min_y = SC_MAX (max_y, y);
      min_z = SC_MAX (max_z, z);

      vertices[3 * v + 0] = x;
      vertices[3 * v + 1] = y;
      vertices[3 * v + 2] = z;

      ++v;

      if (v == num_vertices) {
        axis_scale
          = 1. / SC_MAX (SC_MAX (max_x - min_x, max_y - min_y), max_z - min_z);
      }

      continue;
    } else {
      /* all vertices are read, now read faces */
      retval
        = sscanf (line, "%lu, %lu, %lu, %lu", &v2f, &vid[0], &vid[1], &vid[2]);
      if (retval != 4) {
        P4EST_LERROR ("Premature end of file");
        P4EST_FREE (line);
        return 0;
      }
      if (v2f != 3) {
        P4EST_LERROR ("Unsupported face format in file");
        P4EST_FREE (line);
        return 0;
      }
      faces[f].v0[0] = (vertices[3 * vid[0] + 0] - min_x) * axis_scale;
      faces[f].v0[1] = (vertices[3 * vid[0] + 1] - min_y) * axis_scale;
      faces[f].v0[2] = (vertices[3 * vid[0] + 2] - min_z) * axis_scale;

      faces[f].v1[0] = (vertices[3 * vid[1] + 0] - min_x) * axis_scale;
      faces[f].v1[1] = (vertices[3 * vid[1] + 1] - min_y) * axis_scale;
      faces[f].v1[2] = (vertices[3 * vid[1] + 2] - min_z) * axis_scale;

      faces[f].v2[0] = (vertices[3 * vid[2] + 0] - min_x) * axis_scale;
      faces[f].v2[1] = (vertices[3 * vid[2] + 1] - min_y) * axis_scale;
      faces[f].v2[2] = (vertices[3 * vid[2] + 2] - min_z) * axis_scale;

      ++f;
    }
    P4EST_FREE (line);
  }
  return 1;
}

static int
triangulation_read_off_file (p4est_model_t * m, const char * filename)
{
  int retval = 0;
  FILE *fin = NULL;

  P4EST_GLOBAL_PRODUCTIONF ("Reading connectivity from %s\n", filename);

  fin = fopen (filename, "rb");
  if (fin == NULL) {
    P4EST_LERRORF ("Failed to open %s\n", filename);
    return 0;
  }

  if (!triangulation_read_off_file_stream (m, fin)){
    P4EST_LERRORF ("Failed to read %s: pass 1\n", filename);
    return 0;
  }

  retval = fclose (fin);
  if (retval != 0) {
    P4EST_LERRORF ("Failed to close %s\n", filename);
    return 0;
  }

  return 1;
}

static void
triangulation_setup_model (p4est_model_t ** m, const char * filename)
{
  p4est_model_t  *model = P4EST_ALLOC_ZERO (p4est_model_t, 1);
  model->output_prefix = "triangulation";
  model->conn = p4est_connectivity_new_unitsquare ();
  model->geom = NULL;
  model->intersect = triangulation_intersect_model;
  if (!triangulation_read_off_file (m, filename)) {
    P4EST_LERRORF ("Failed to read a valid model from %s\n", filename);
    sc_abort ();
  }
}
/* #endif */

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