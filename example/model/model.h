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

#ifndef P4EST_MODEL_H
#define P4EST_MODEL_H

#include <p4est.h>
#include <p4est_geometry.h>

/** Used to free private model primitives. */
typedef void        (*p4est_destroy_primitives_t) (void *primitives);

/** Check intersection of a rectangle with an object. */
typedef int         (*p4est_intersect_t) (p4est_topidx_t which_tree,
                                          const double coord[4],
                                          void *primitive);

/** General, application specific model data.
 * Creating the model as well as design and managing of primitives is
 * completely up to a user.
 */
typedef struct p4est_model
{
  const char         *output_prefix;
  size_t              num_prim;   /**< number of primitives in the model */
  p4est_connectivity_t *conn;     /**< connectivity the model exists in */
  p4est_geometry_t   *geom;       /**< deformation parameters for output */
  void               *primitives; /**< ptr to stored model's primitives */

  /** When not NULL, free whatever is stored in model::primitives. */
  p4est_destroy_primitives_t destroy_primitives;

  /** Intersect a given rectangle with a model primitive. */
  p4est_intersect_t intersect;

  /** Private geometry data. */
  /*p4est_geometry_t    sgeom;*/
}
p4est_model_t;

/** Wrapping for model's functions */
/** Destroy model */
void                p4est_model_destroy (p4est_model_t * model);

/** Implementation of p4est's callbacks w.r.t. the model structure */
/** Callback function to query the match of a "point" with a quadrant.
 * Its prototype completely corresponds to p4est_search_local_t
 * type from p4est_search.h.
 */
int                 p4est_model_intersect (p4est_t * p4est,
                                           p4est_topidx_t which_tree,
                                           p4est_quadrant_t * quadrant,
                                           p4est_locidx_t local_num,
                                           void *point);

/** Callback function prototype to decide for refinement.
 * Its prototype completely corresponds to p4est_refine_t type from p4est.h.
*/
int                 p4est_model_refine (p4est_t * p4est,
                                        p4est_topidx_t which_tree,
                                        p4est_quadrant_t *quadrant);

/** Callback function prototype to initialize the quadrant's user data.
 * Its prototype completely corresponds to p4est_init_t type from p4est.h.
*/
void                p4est_model_quad_init (p4est_t * p4est,
                                           p4est_topidx_t which_tree,
                                           p4est_quadrant_t *quadrant);

#endif /* P4EST_MODEL_H */
