
# This file is part of p4est.
# Makefile.am in example/gmt
# included non-recursively from toplevel directory

if P4EST_ENABLE_BUILD_2D
bin_PROGRAMS += example/model/p4est_refsearch

example_model_p4est_refsearch_SOURCES = \
        example/model/model.c \
        example/model/refsearch2.c

LINT_CSOURCES += \
        $(example_model_p4est_refsearch_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
endif
