SetFactory("OpenCASCADE");
lc = 0.1;

// Define corner points of L-shaped domain
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.5, 0, lc};
Point(4) = {0.5, 0.5, 0, lc};
Point(5) = {0.5, 1, 0, lc};
Point(6) = {0, 1, 0, lc};

// Define lines (OpenCASCADE-compatible)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// Create a wire (closed loop) for the L-shape
Wire(1) = {1, 2, 3, 4, 5, 6};

// Create the L-shaped surface (OpenCASCADE face)
Surface(1) = {1};

// Create a circular disk (to subtract)
Disk(2) = {0.25, 0.25, 0, 0.2, 0.2};  // center (x, y, z), radius_x, radius_y

// Subtract the disk from the L-shape
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }

// Assign physical group to final surface
Physical Surface("L_with_hole") = {1};

