Point(1) = {0, 0, 0, 1};
Extrude {1, 0, 0} {
  Point{1};Layers{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{1};Layers{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{2};Layers{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{5};Layers{4};
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Line{8};Layers{4};
}
