// Inputs
	lc=1;
	lx=1.0000000000000000000000000000000;
	ly=1.0000000000000000000000000000000;
	x0=0.1;
	x1=0.9;
	
	y0=0.1;
	y1=0.9; 	

        // Geometry
	Point(1) = {0, 0, 0, lc};
	Point(2) = {x0*lx, 0, 0, lc};
	Point(3) = {x1*lx, 0, 0, lc};
	Point(4) = {lx, 0, 0, lc};

	Point(5) = {lx, y0*ly, 0, lc};
	Point(6) = {lx, y1*ly, 0, lc};

	Point(7) = {lx, ly, 0, lc};
	Point(8) = {x1*lx, ly, 0, lc};
	Point(9) = {x0*lx, ly, 0, lc};
	Point(10) = {0, ly, 0, lc};
	Point(11) = {0, y1*ly, 0, lc};
	Point(12) = {0, y0*ly, 0, lc};

	Line(1) = {1, 2};				
	Line(2) = {2, 3};				
	Line(3) = {3, 4};				
	Line(4) = {4, 5};
	Line(5) = {5, 6};				
	Line(6) = {6, 7};				
	Line(7) = {7, 8};				
	Line(8) = {8, 9};
	Line(9) = {9, 10};				
	Line(10) = {10, 11};				
	Line(11) = {11, 12};				
	Line(12) = {12, 1};
	
	nxb=10;
	nxm=20;
	nyb=10;
	nym=20;
	Transfinite Line(1)=	nxb;
	Transfinite Line(3)=	nxb;
	Transfinite Line(7)=	nxb;
	Transfinite Line(9)=	nxb;

	Transfinite Line(2)=	nxm;
	Transfinite Line(8)=	nxm;
	
	Transfinite Line(4)=	nyb;
	Transfinite Line(6)=	nyb;
	Transfinite Line(10)=	nyb;
	Transfinite Line(12)=	nyb;

	Transfinite Line(5)=	nym;
	Transfinite Line(11)=	nym;
		
	Line Loop(5) = {1, 2, 3, 4,5,6,7,8,9,10,11,12};
	//Line Loop(7) = {4,1,2};  	
	Plane Surface(6) = {5};
 
        //Transfinite surface:
	Transfinite Surface {6};
	Recombine Surface {6};

	Physical Line("driven_lid",1)={3};
	Physical Line("wall",2)={4,1,2};
	Physical Surface("flow_field",200)={6};
 

