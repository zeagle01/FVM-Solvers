#ifndef LAPLACE_SOLVER_H_H
#define LAPLACE_SOLVER_H_H



#include "solver.h"
#include "Laplace.h"
#include "dataStructure.h"
#include "configReader.h"
#include "serviceLocator.h"

class DiffusionSolver :public Solver, public member_t<DiffusionSolver>{
public:

	I_Laplace* diffusion;
	double gama;

	virtual void config(){
		Solver::config();
		gama = configReader->readDouble("gama");
		diffusion=configReader->readLaplace("Laplace_operator");
	}

	virtual void solve(){
		config();
		
		CSR eq = diffusion->apply(gama, phi[0], mesh);
		//les[0]->solve(eq, phi[0], mesh);
		les[0]->solve(eq, phi[0], mesh);
		phi[0].assignBoundary(mesh);
		printers[0]->print("out.dat",mesh, phi);
		
	}
	
};

/*
class UnorthogonalCorrectionDiffusionSolver :public Solver{
public:

	I_Laplace* diffusion;
	double gama;
	I_GaussCellGradient* grad;
	vector<I_FaceReconstruct*> faceReconstruct;
	UnorthogonalCorrectionDiffusionSolver(Case* cas) :Solver(cas) {

	}

	virtual void solve(){
		cas->config();
		gama = cas->configReader->readDouble("gama");
		diffusion = cas->configReader->readLaplace("Laplace_operator");
		grad = cas->configReader->readCellGradient("cell_gradient");

		faceReconstruct = cas->configReader->readFaceReconstructor("face_average");
		for (int it = 0; it < 100; it++){
			vector<double> phi_f = faceReconstruct[0]->apply(cas->phi[0], cas->mesh);
			vector<CellField> grad_c = grad->apply(phi_f, cas->mesh);
			vector<double> grad_f_x = faceReconstruct[1]->apply(grad_c[0], cas->mesh);
			vector<double> grad_f_y = faceReconstruct[1]->apply(grad_c[1], cas->mesh);

			vector<double> grad_f = grad_f_y;
			grad_f.insert(grad_f.begin(), grad_f_x.begin(), grad_f_x.end());

			diffusion->grad_f = grad_f;
			CSR eq = diffusion->apply(gama, cas->phi[0], cas->mesh);
			cas->les[0]->solve(eq, cas->phi[0], cas->mesh);
			cas->phi[0].assignBoundary(cas->mesh);
			cas->printers[0]->print("out.dat", cas->mesh, cas->phi);
		}


		//std::ofstream fout("a.dat");
		//copy(phi[0]->inner.cbegin(), phi[0]->inner.cend(), ostream_iterator<T>(fout, "\n") );
		//for_each(phi[0]->inner.begin(), phi[0]->inner.end(), print<T>());
	}
};


*/



#endif