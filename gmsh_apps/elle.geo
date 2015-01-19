// Gmsh project created on Mon Jan 19 09:49:17 2015
Point(1) = {0, 0, 0, .2};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {0, .5, 0} {
  Point{2};
}
Extrude {-0.5, 0, 0} {
  Point{3};
}
Extrude {0, 0.5, 0} {
  Point{4};
}
Extrude {-0.5, 0, 0} {
  Point{5};
}
Extrude {0, -1, 0} {
  Point{6};
}
Line Loop(7) = {6, 1, 2, 3, 4, 5};
Plane Surface(8) = {7};
