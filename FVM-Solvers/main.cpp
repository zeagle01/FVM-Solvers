
#include "solver.h"
#include "solverRegister.h"

#include <boost/timer.hpp>
#include <iostream>
using namespace std;
int main(){

	cout << "input the CONFIG FILE(.json),press enter to finish the input" << endl;
	cout << "or just drag the CONFIG FILE(.json) into the console" << endl;
	string s;
	getline(cin,s);
	//s = "square_cylinder.json";
	boost::timer t;  //定义一个计时类，开始计时
	SolverRegister solverRegister(s);
	Solver* solver = solverRegister.readSolver("solver");
	cout << "compute begins" << endl;
	solver->solve();

	std::cout << "computation time：" << t.elapsed() << std::endl;
	system("pause");
	//no ff merge;
	return 0;
}