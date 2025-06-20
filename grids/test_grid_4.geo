SetFactory("OpenCASCADE");

// Outer polygon
Point(1) = {0, 1, 0, 1.0};
Point(2) = {1, 1, 0, 1.0};
Point(3) = {2, 0, 0, 1.0};
Point(4) = {1, -1, 0, 1.0};
Point(5) = {0, -1, 0, 1.0};
Point(6) = {-1, 0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {5, 4};
Line(5) = {6, 5};
Line(6) = {1, 6};

Curve Loop(10) = {1, 2, -3, 4, 5, 6};

// Hole (circular arc)
Point(7) = {0.5, 0.75, 0, 1.0};       // Center
Point(8) = {0.75, 0.75, 0, 1.0};      // Start of arc
Point(9) = {0.25, 0.75, 0, 1.0};      // End of arc

// Use Circle with center point for OpenCASCADE (creates trimmed arc)
Circle(20) = {8, 7, 9};     // Arc from 8 -> 9 with center 7
Circle(21) = {9, 7, 8};     // Close the arc (other direction to form full loop)

Curve Loop(22) = {20, 21};  // Loop representing full circle (can be a disk hole)

// Create face with hole
Plane Surface(11) = {10, 22};

// Optionally mark boundaries
Physical Surface("Domain") = {11};
Physical Curve("Outer") = {1,2,3,4,5,6};
Physical Curve("ArcHole") = {20,21};

//+
Show "*";
//+
Show "*";
