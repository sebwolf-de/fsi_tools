cl1=1;
ly=2;
lx=4;
Point(1)={0,0,0,cl1};
Extrude {1,0,0}
{
Point{1};
Layers{lx};
}
Extrude{0,1,0}
{
Line{1};
Layers{ly};
}
