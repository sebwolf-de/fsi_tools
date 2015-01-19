// Gmsh project created on Wed Jan 14 10:22:05 2015
Point(1) = {0, 0, 0, 0.5};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{1};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{2};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{5};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{8};
}
