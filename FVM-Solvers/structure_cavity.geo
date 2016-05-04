lx=1.0;
ly=1.0;
<<<<<<< HEAD
nx=51;
ny=51;
=======
nx=101;
ny=101;
>>>>>>> dev

bump_x=1.0;
bump_y=1.0;


Point(1)=(0,0,0);
Point(2)=(lx,0,0);
Point(3)=(lx,ly,0);
Point(4)=(0,ly,0);


Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};



Transfinite Line{1,-3}= nx Using Bump bump_x;
Transfinite Line{2,-4}= ny Using Bump bump_y;

	Line Loop(5) = {1, 2, 3, 4};
		
	Plane Surface(5) = {5};
	
 
        //Transfinite surface:
	Transfinite Surface {5};
	Recombine Surface {5};


Physical Line("driven_lid",1)={3};
	Physical Line("wall",2)={4,1,2};
	Physical Surface("flow_field",200)={5};

