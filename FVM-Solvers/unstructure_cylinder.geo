lc = 0.1;
lx=2.0;
ly=2.0;
Point(1) = {-lx, -ly, 0.0, lc};
Point(2) = { lx, -ly, 0.0, lc};
Point(3) = { lx,  ly, 0.0, lc};
Point(4) = {-lx,  ly, 0.0, lc};

r=0.8;
Point(5) = {0.0, 0.0, 0.0, lc};
Point(6) = {r, 0.0, 0.0, lc};
Point(7) = {0.0, r, 0.0, lc};
Point(8) = {-r, 0.0, 0.0, lc};
Point(9) = {0.0, -r, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
nx=41;ny=41;
Transfinite Line{1} = nx;
Transfinite Line{2} = ny;
Transfinite Line{3} = nx;
Transfinite Line{4} = ny;

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

n_cylinder=41;
Transfinite Line{5} = n_cylinder;
Transfinite Line{6} = n_cylinder;
Transfinite Line{7} = n_cylinder;
Transfinite Line{8} = n_cylinder;

Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5,6,7,8};
Plane Surface(13) = {11,12};


Physical Line("outSide",1)={1,2,3,4};
Physical Line("cylinder",2)={5,6,7,8};
Physical Surface("flow_field",200)={13};

