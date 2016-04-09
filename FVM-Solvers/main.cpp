
#include "solver.h"
#include "solverRegister.h"

#include <boost/timer.hpp>
int main(){

	boost::timer t;  //定义一个计时类，开始计时
	string s = "pressure_correction_solver.json";
	SolverRegister solverRegister(s);
	Solver* solver = solverRegister.readSolver("solver");
	solver->solve();

	std::cout << "运行时间：" << t.elapsed() << std::endl;
	system("pause");
	return 0;
}