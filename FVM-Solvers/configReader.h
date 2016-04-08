#ifndef CONFIGREADER_H_H
#define CONFIGREADER_H_H

#include <boost/foreach.hpp>
#include <boost\property_tree\json_parser.hpp>
#include <string>
#include "serviceLocator.h"
#include "mesh.h"
#include "laplace.h"
#include "linearEquationSolver.h"
#include "printer.h"
using namespace std;



class A_ConfigReader{
	string configFile;
	boost::property_tree::ptree pt;
	servicelocator_t services;
public:
	A_ConfigReader(string configFile){
		this->configFile = configFile;
		boost::property_tree::json_parser::read_json(configFile, pt);
		buildService();
	}
	//register the class
	void buildService(){
		services.register_class<GmeshTriangleMesh>("Gmesh_triangle_mesh");
		services.register_class<VersatileGmeshMesh>("versatile_Gmesh_mesh");
		services.register_class<Laplace>("Laplace");
		services.register_class<GS_Solver>("GS_iteration");
		services.register_class<UnstructurePrinter>("tecplot_unstructure");
	}

	int readInt(string s){
		return pt.get_child(s).get_value<int>();
	}
	double readDouble(string s){
		return pt.get_child(s).get_value<double>();
	}

	string readString(string s){
		return pt.get_child(s).get_value<string>();
	}



	Mesh* readMeshClass(string s){
		Mesh* r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_child("class").get_value<string>();
		string meshFile = child.get_child("mesh_file").get_value<string>();
		r = services.get_single_instance<Mesh>(type);
		r->meshFile = meshFile;
		return r;
	}

	template<class T> vector<T> read_1D_array(string s){

		vector<T> r;
		boost::property_tree::ptree child_array = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child_array) {
			r.push_back(vt.second.get_value<T>());
		}
		return r;
	}

	I_Laplace* readLaplace(string s){
		I_Laplace* r;
		string type = pt.get_child(s).get_value<string>();
		r = services.get_single_instance<I_Laplace>(type);
		return r;
	}


	vector<I_LinearEquationSolver*> readLinearEquationSolver(string s){
		vector<I_LinearEquationSolver*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			r.push_back(services.get_single_instance<I_LinearEquationSolver>(type));
			r.back()->converge_threhold = vt.second.get_child("converge_threhold").get_value<double>();
			r.back()->max_step = vt.second.get_child("max_step").get_value<int>();
			r.back()->check_step = vt.second.get_child("check_step").get_value<int >();
		}
		return r;
	}

	vector<A_Printer*> readPrinter(string s){
		vector<A_Printer*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			string zone_type = vt.second.get_child("zone_type").get_value<string>();

			r.push_back(services.get_single_instance<A_Printer>(type));
			r.front()->zoneType = zone_type;
		}
		return r;
	}

};

#endif