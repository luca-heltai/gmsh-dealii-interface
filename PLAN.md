# Plan of action

Serial version (currently working):

- [X] create a simple grid with gmsh (.geo file)
- [X] make a (single) `.gmsh` file with gmsh
- [x] read the `.gmsh` file with gmsh using the existing gmsh API in dealii

Parallel version (not working yet):

- [ ] create as many `.msh` files as there are processes
- [ ] read the `.msh` files with gmsh using the existing gmsh API in dealii separately in each process
- [ ] merge the meshes in each process using the existing gmsh API
- [ ] write a code that generates a parallel::fullydistributed::Triangulation from the merged mesh
