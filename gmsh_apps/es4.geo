cl1=.01;
cl2=.4;

Point(1)={0,0,0,cl1};
Point(2)={1,0,0,cl1};
Point(3)={1,1,0,cl1};
Point(4)={0,1,0,cl1};
Point(5)={0.5,0,0,cl2};
Point(6)={1,.5,0,cl2};
Point(7)={.5,1,0,cl2};
Point(8)={0,.5,0,cl2};

Line(1)={1,5};
Line(2)={5,2};
Line(3)={2,6};
Line(4)={6,3};
Line(5)={3,7};
Line(6)={7,4};
Line(7)={4,8};
Line(8)={8,1}; 

Line Loop(5)={1,2,3,4,5,6,7,8};

Plane Surface(6)={5};

