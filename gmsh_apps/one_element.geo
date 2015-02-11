// Gmsh project created on Wed Feb 11 15:05:48 2015
Point(1) = {0, 0, 0, 2.0};
Point(2) = {1, 0, 0, 2.0};
Point(3) = {0, 1, 0, 2.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};
