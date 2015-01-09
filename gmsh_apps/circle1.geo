Point(1) = {0, 0, 0, .2};
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
