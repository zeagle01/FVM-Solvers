#ifndef SOLVER_REGISTER_H_H
#define SOLVER_REGISTER_H_H

#include <boost/foreach.hpp>
#include <boost\property_tree\json_parser.hpp>
#include <string>
#include "serviceLocator.h"
#include "solver.h"
#include "laplaceSolver.h"
using namespace std;


class SolverRegister{
	string configFile;
	boost::property_tree::ptree pt;
	servicelocator_t services;
public:
	SolverRegister(string configFile){
		this->configFile = configFile;
		boost::property_tree::json_parser::read_json(configFile, pt);
		buildService();
	}
	//register the class
	void buildService(){
		services.register_class<DiffusionSolver>("laplace_solver");
	}

	Solver* readSolver(string s){
		Solver* r;
		string type = pt.get_child(s).get_value<string>();
		r = services.get_single_instance<Solver>(type);
		r->configReader = new A_ConfigReader(configFile);
		return r;
	}

};

#endif