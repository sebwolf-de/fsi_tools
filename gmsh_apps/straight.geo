lc = 0.1;

Point(1) = {0,    0,    0, lc};
Point(2) = {2.5,  0,    0, lc};
Point(3) = {2.5,  0.41, 0, lc};
Point(4) = {0  ,  0.41, 0, lc};

lc = 0.05;

Point(5) = {0.2,  0.2,  0, lc};

Point(6) = {0.25, 0.2,  0, lc};
Point(7) = {0.2,  0.25, 0, lc};
Point(8) = {0.15, 0.2,  0, lc};
Point(9) = {0.2,  0.15, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Line Loop(9) = {1,2,3,4}; // exterior loop
Line Loop(10) = {5,6,7,8};

Plane Surface(1) = {9,10}; // exterior surface (with a whole)
