lc1 = 2;
lc2 = 1;
lc3 = 0.5;
lc4 = 0.25;

f = 0.125;

Point(1) = {0,  0,    0, f*lc3};
Point(2) = {5,  0,    0, f*lc2};
Point(3) = {10, 0,    0, f*lc1};
Point(4) = {10, 1.64, 0, f*lc1};
Point(5) = {5,  1.64, 0, f*lc2};
Point(6) = {0,  1.64, 0, f*lc3};

//Point(5) = {0.8, 0.8,  0};

Point(10) = {0.6, 0.6, 0, f*lc4};
Point(11) = {1,   0.6, 0, f*lc4};
Point(12) = {1,   1,   0, f*lc4};
Point(13) = {0.6, 1,   0, f*lc4};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};
//Circle(5) = {6,5,7};
//Circle(6) = {7,5,8};
//Circle(7) = {8,5,9};
//Circle(8) = {9,5,6};
Line Loop(9) = {1, 2, 3, 4, 5, 6}; // exterior loop
Line Loop(10) = {10, 11, 12, 13};  // interior loop
//Plane Surface(1) = {9}; // interior surface
Plane Surface(2) = {9,10}; // exterior surface (with a whole)
