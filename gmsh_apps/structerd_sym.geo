cl1 = 0.1;
lyx = 3;
Point(1) = {0, 0.5, 0, cl1};
Extrude {1, 0, 0}
{
  Point{1};Layers{lyx};
}
Extrude {0, .5, 0}
{
  Line{1};Layers{lyx};
}
Extrude {0, -.5, 0}
{
  Line{1};Layers{lyx};
}
