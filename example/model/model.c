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

#include "model.h"

static const double irootlen = 1. / (double) P4EST_ROOT_LEN;

void
p4est_model_destroy (p4est_model_t * model)
{
  if (model->destroy_primitives != NULL) {
    model->destroy_primitives(model->primitives);
  }
}

int
p4est_model_intersect (p4est_t * p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t * quadrant, p4est_locidx_t local_num,
                       void *point)
{
  int                 result;
  double              coord[6];
  size_t              p;
  p4est_qcoord_t      qh;
  p4est_model_t  *model = (p4est_model_t *) p4est->user_pointer;

  /* sanity checks */
  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (model->intersect != NULL);

  /* retrieve object index and model */
  P4EST_ASSERT (point != NULL);
  p = *(size_t *) point;
  P4EST_ASSERT (p < model->num_prim);

  /* provide rectangle coordinates */
  qh = P4EST_QUADRANT_LEN (quadrant->level);
  coord[0] = irootlen * quadrant->x;
  coord[1] = irootlen * quadrant->y;
#ifdef P4EST_TO_P8EST
  coord[2] = irootlen * quadrant->z;
#endif
  coord[3] = irootlen * (quadrant->x + qh);
  coord[4] = irootlen * (quadrant->y + qh);
#ifdef P4EST_TO_P8EST
  coord[5] = irootlen * (quadrant->z + qh);
#endif

  /* execute intersection test */
  if ((result = model->intersect (which_tree, coord, model, point)) &&
      local_num >= 0) {
    /* set refinement indicator for a leaf quadrant */
    quadrant->p.user_int = 1;
  }
  return result;
}

int
p4est_model_refine (p4est_t * p4est, p4est_topidx_t which_tree,
                    p4est_quadrant_t *quadrant)
{
  return quadrant->p.user_int;
}

void
p4est_model_quad_init (p4est_t * p4est, p4est_topidx_t which_tree,
                       p4est_quadrant_t *quadrant)
{
  quadrant->p.user_int = 0;
}
