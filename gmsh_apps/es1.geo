cl1 = 2;
lyx = 5;
Point(1) = {0, 0.0, 0, cl1};
Extrude {1, 0, 0} 
{
  Point{1};
  Layers{lyx};
}
Extrude {0, 1, 0} 
{
  Line{1};
  Layers{lyx};
}
