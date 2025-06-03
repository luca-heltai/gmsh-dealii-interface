// Gmsh project created on Tue Jun  3 12:16:45 2025
SetFactory("OpenCASCADE");

c = .3; // Characteristic length

//+
Point(1) = {-1, 1, 0, c};
//+
Point(2) = {0, 0, 0, c};
//+
Line(1) = {1, 2};
//+
Extrude {{0, 1, 0}, {1, -1, 0}, Pi/2} {
  Curve{1}; 
}
//+
Physical Surface(5) = {1};
