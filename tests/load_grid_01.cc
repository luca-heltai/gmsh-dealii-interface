// ---------------------------------------------------------------------
//
// Copyright (C) 2024 by Luca Heltai
//
// This file is part of the bare-dealii-app application, based on the
// deal.II library.
//
// The bare-dealii-app application is free software; you can use it,
// redistribute it, and/or modify it under the terms of the Apache-2.0 License
// WITH LLVM-exception as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md
// at the top level of the bare-dealii-app distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "gmsh_api.h"
#include "tests.h"

auto
main() -> int
{
  initlog(1);

  // Create a 2D triangulation
  dealii::Triangulation<2, 3> triangulation;

  // Open and read the mesh file
  std::string mesh_file = SOURCE_DIR "/../grids/test_grid_1.msh";

  GMSH::read_msh(triangulation, mesh_file);

  // Output the number of cells and vertices
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  deallog << "Number of vertices: " << triangulation.n_vertices() << std::endl;

  return 0;
}
