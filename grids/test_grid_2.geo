//+
SetFactory("OpenCASCADE");
//+
lc = 0.1;
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {1, 1, 0, lc};
//+
Point(4) = {0, 1, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(6) = {5};
//+
Physical Point("Boundary Points") = {1, 2, 3, 4}; // Physical group for points
Physical Line("Boundary Lines") = {1, 2, 3, 4};   // Physical group for lines
Physical Surface("Rectangle Surface") = {6};      // Physical group for the surface
//+
Line(5) = {4, 2};
//+
Point(5) = {-0, 0.5, 0, lc};
//+
Point(6) = {1, 0.5, 0, lc};
//+
Point(7) = {0, 0.5, 0, lc};
//+
Line(6) = {6, 4};
