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
#include "faceReconstruct.h"
#include "convection.h"
#include "temporalTerm.h"
#include "gradient.h"
#include "sourceTerm.h"
#include "explicitDivergence.h"
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
		services.register_class<UpwindConvection>("1st_upwind");
		services.register_class<DefferCenterConvection>("DC_center");
		services.register_class<ExplicitCenter>("explicit_center");
		services.register_class<ImplicitCenter>("implicit_center");
		services.register_class<ImplictEuler>("1st_Euler");
		services.register_class<GS_Solver>("GS_iteration");
		services.register_class<UnstructurePrinter>("tecplot_unstructure");
		services.register_class<FaceReconstruct_BEB>("linear_average_beb");
		services.register_class<LinearFaceReconstruct_BEO>("linear_average_beo");
		services.register_class<GaussCellGradient>("GS_gradient");
		services.register_class<PlainFaceGradient>("plain_face_gradient");
		services.register_class<SourceTerm>("field_source_term");
		services.register_class<RhieChowInterpolation>("basic_Rhie_Chow");
		services.register_class<ExplicitDivergence>("basic_explicit_divergence");
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
	bool readBool(string s){	
		return pt.get_child(s).get_value<bool>();
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

	
		I_Convection* readConvection(string s){
			I_Convection* r;
			string type = pt.get_child(s).get_value<string>();
			r = services.get_single_instance<I_Convection>(type);
			return r;
		}


		I_Temporal* readTemporal(string s){
			I_Temporal* r;
			string type = pt.get_child(s).get_value<string>();
			r = services.get_single_instance<I_Temporal>(type);
			r->dt = this->readDouble("dt");
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



	I_GaussCellGradient* readCellGradient(string s){
		I_GaussCellGradient *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<I_GaussCellGradient>(type);
		return r;
	}


	I_FaceGradient* readFaceGradient(string s){
		I_FaceGradient *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<I_FaceGradient>(type);
		return r;
	}


	vector<I_FaceReconstruct*> readFaceReconstructor(string s){
		vector<I_FaceReconstruct*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			r.push_back(services.get_single_instance<I_FaceReconstruct>(type));
		}
		return r;
	}

	I_SourceTerm* readSourceTerm(string s){
		I_SourceTerm *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<I_SourceTerm>(type);
		return r;
	}
	B_RhieChow* readRhieChowInterpolation(string s){
		B_RhieChow *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<B_RhieChow>(type);
		r->faceReconstruct = readFaceReconstructor("face_average");
		return r;
	}
	I_ExplicitDivergence* readDivergence(string s){
		I_ExplicitDivergence *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<I_ExplicitDivergence>(type);
		return r;
	}
};

#endif