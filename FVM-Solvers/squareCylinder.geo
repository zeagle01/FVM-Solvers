// Inputs
	lc=0.1;
	lx=1.0;
	ly=1.0;
 	
	lxi=0.2;
	lyi=0.3;
      xci=0.5;
	yci=0.5;  

	Point(1) = {0, 0, 0, lc};
	Point(2) = {lx, 0, 0, lc};
	Point(3) = {lx, ly, 0, lc};
	Point(4) = {0, ly, 0, lc};
        
	Point(5) = {xci*lx, yci*ly, 0, lc};
	Point(6) = {(xci+lxi)*lx, yci*ly, 0, lc};
	Point(7) = {(xci+lxi)*lx, (yci+lyi)*ly, 0, lc};
	Point(8) = {xci*lx, (yci+lyi)*ly, 0, lc};



	Line(1) = {1, 2};				// bottom line
	Line(2) = {2, 3};				// right line
	Line(3) = {3, 4};				// top line
	Line(4) = {4, 1};				// left line
	Line(5) = {5, 6};				// bottom line
	Line(6) = {6, 7};				// right line
	Line(7) = {7, 8};				// top line
	Line(8) = {8, 5};				// left line


	Line Loop(5) = {1, 2, 3, 4};
	Line Loop(6) = {5, 6, 7, 8};
	//Line Loop(7) = {4,1,2};  	
	Plane Surface(5) = {5};
	//Plane Surface(6) = {6};
 
        //Transfinite surface:
	Transfinite Surface {5};
	Recombine Surface {5};

	//Transfinite Surface {6};
	//Recombine Surface {6};

	//Physical Line("driven_lid",1)={3};
	//Physical Line("wall",2)={4,1,2};
	//Physical Surface("flow_field",200)={6};
