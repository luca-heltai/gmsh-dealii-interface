SetFactory("OpenCASCADE");
lc = 0.1;

// Define points
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Define lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Define physical groups
Physical Surface("Rectangle Surface") = {6}; // Only the surface is relevant for partitioning
