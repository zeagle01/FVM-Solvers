// Inputs
	squareSide = 1; //m
	//meshThickness = squareSide / 10; 
	gridsize = squareSide/50 ;
 
        // Geometry
	Point(1) = {-squareSide/2, -squareSide/2, 0, gridsize};
	Point(2) = {squareSide/2, -squareSide/2, 0, gridsize};
	Point(3) = {squareSide/2, squareSide/2, 0, gridsize};
	Point(4) = {-squareSide/2, squareSide/2, 0, gridsize};
	Line(1) = {1, 2};				// bottom line
	Line(2) = {2, 3};				// right line
	Line(3) = {3, 4};				// top line
	Line(4) = {4, 1};				// left line
	Line Loop(5) = {1, 2, 3, 4};
	//Line Loop(7) = {4,1,2};  	
	Plane Surface(6) = {5};
 
        //Transfinite surface:
	//Transfinite Surface {6};
	//Recombine Surface {6};

	Physical Line("driven_lid",1)={3};
	Physical Line("wall",2)={4,1,2};
	Physical Surface("flow_field",200)={6};
 

