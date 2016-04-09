// Inputs
	nx=50;
	ny=50;
 	
	nh=10;
	nw=20;
	
	lx=1.0;
	ly=1.0;
	h=ly*(nh/ny);
	w=lx*(nw/nx);
	Point(1) = {0, h, 0, 1};
	Point(2) = {w, h, 0, 1};
	Point(3) = {w, 0, 0, 1};
	Point(4) = {lx, 0, 0, 1};
	Point(5) = {lx, ly, 0, 1};
	Point(6) = {0, ly, 0, 1};


	Line(1) = {1, 2};				
	Line(2) = {2, 3};				
	Line(3) = {3, 4};				
	Line(4) = {4, 5};				
	Line(5) = {5, 6};				
	Line(6) = {6, 1};				

	Transfinite Line{1} = nw;
	Transfinite Line{2} = nh;	
	Transfinite Line{3} = nx-nw;
	Transfinite Line{4} = ny;	
	Transfinite Line{5} = nx;
	Transfinite Line{6} = ny-nh;		
	
	Line Loop(5) = {1, 2, 3, 4,5,6};
		
	Plane Surface(5) = {5};
	
 
        //Transfinite surface:
	Transfinite Surface {5};
	Recombine Surface {5};


	Physical Line("left",1)={6};
	Physical Line("right",2)={4};
	Physical Line("wall",3)={1,2,3,5};

	Physical Surface("flow_field",200)={5};
