#ifndef PRESSURE_CORRECTION_SOLVER_H_H
#define PRESSURE_CORRECTION_SOLVER_H_H

#include "faceReconstruct.h"
#include "convection.h"
#include "laplace.h"
#include "dataStructure.h"
#include "solver.h"
#include "temporalTerm.h"
#include "gradient.h"
#include "sourceTerm.h"
#include "explicitDivergence.h"

class PressureCorrectionSolver :public Solver, public member_t<PressureCorrectionSolver>{
public:
	I_Laplace* diffusion;
	I_Convection* convection;
	I_Temporal* temporal;

	I_GaussCellGradient* cGrad;
	I_FaceGradient* fGrad;
	I_ExplicitDivergence* div;
	vector<I_FaceReconstruct*> faceReconstruct;
	B_RhieChow* rc;
	I_SourceTerm* source;

	bool relativePressure;

	vector<CellField*> V;

	double Re;
	double gama;
	double dt;
	int step;
	int maxStep;
	int checkStep;
	double  error;
	double convergenceThrehold;

	virtual void config(){
		Solver::config();
		Re = configReader->readDouble("Re");
		diffusion = configReader->readLaplace("Laplace_operator");
		convection = configReader->readConvection("convection_operator");
		faceReconstruct = configReader->readFaceReconstructor("face_average");
		temporal = configReader->readTemporal("temporal_scheme");
		cGrad = configReader->readCellGradient("cell_gradient");
		fGrad = configReader->readFaceGradient("face_gradient");
		source = configReader->readSourceTerm("source_term_operator");
		rc = configReader->readRhieChowInterpolation("Rhie_Chow");
		div = configReader->readDivergence("explicit_divergence");
		relativePressure = configReader->readBool("relative_pressure");
		dt = temporal->dt;


		V.push_back(&(phi[0]));
		V.push_back(&(phi[1]));
		
		gama = 1.0 / Re;

		checkStep = configReader->readInt("check_step");
		maxStep = configReader->readInt("max_step");
		convergenceThrehold = configReader->readDouble("converge_threhold");

	}
	virtual void solve(){
		config();
		vector<CSR> V_diffusion = diffusion->apply(gama, V, mesh);
		CSR  p_eq = diffusion->apply(1.0, phi[2], mesh);
		CellField p_prime = phi[2];
		int step = 0;
		vector<CellField> phi0(3, CellField());
		do{
			if (step%checkStep == 0){
				phi0[0] = phi[0];
				phi0[1] = phi[1];
				phi0[2] = phi[2];
			}
			vector<double> V_f = faceReconstruct[0]->apply(V, mesh);




			vector<CSR> V_convection = convection->apply(V_f, V, mesh);
			vector<CSR> V_timeTerm = temporal->apply(V, mesh);
			vector<CSR> V_momentum = CSR::minus(V_convection, V_diffusion);
			V_momentum = CSR::plus(V_momentum, V_timeTerm);
			

			vector<double> p_f = faceReconstruct[0]->apply(phi[2],mesh);
			vector<double> p_cellGrad = cGrad->apply(p_f, mesh);
			vector<CSR> sourceTerm = source->apply(V, mesh, p_cellGrad);

			V_momentum = CSR::plus(V_momentum, sourceTerm);

			//les[0]->iterate(V_momentum[0], phi[0], mesh, 1);
			//les[0]->iterate(V_momentum[1], phi[1], mesh, 1);
			les[0]->solve(V_momentum[0], phi[0], mesh);
			les[0]->solve(V_momentum[1], phi[1], mesh);
			phi[0].assignBoundary(mesh);
			phi[1].assignBoundary(mesh);

		
			vector<double> fgrad = fGrad->apply(phi[2], mesh);

			vector<double> D_f = rc->calculateD(V_momentum, mesh);

			vector<CellField> p_cellGradField = CellField::wrapInner(p_cellGrad,2,mesh);
			vector<double> v_Rc = rc->reconstructVelocity(D_f, p_cellGradField, fgrad, V, mesh);
			vector<double> mdot = div->apply(v_Rc, mesh);

			for (int c = 0; c<mesh->cellNum; c++){
				p_eq.b[c] = mdot[c] * mesh->C_v[c] / dt;
			}

			//les[1]->iterate(p_eq, p_prime, mesh, 1);
			les[1]->solve(p_eq, p_prime, mesh);
			
			if (relativePressure == true){
				setRelativePressure(p_prime);
			}
			
			p_prime.assignBoundary(mesh);

			vector<double> p_prime_f = faceReconstruct[0]->apply(p_prime, mesh);
			vector<double> p_prime_cgrad = cGrad->apply(p_prime_f, mesh);
			for (int c = 0; c < mesh->cellNum; c++) {
				phi[0].inner[c] = phi[0].inner[c] - p_prime_cgrad[c] * dt;
				phi[1].inner[c] = phi[1].inner[c] - p_prime_cgrad[c+mesh->cellNum] * dt;
				phi[2].inner[c] += p_prime.inner[c];
			}
			phi[0].assignBoundary(mesh);
			phi[1].assignBoundary(mesh);
			phi[2].assignBoundary(mesh);

			if (step%checkStep == 0){
				error = VectorMath<double>::rootOfSquareSum(phi[0].inner, phi0[0].inner);
				error += VectorMath<double>::rootOfSquareSum(phi[1].inner, phi0[1].inner);
				error += VectorMath<double>::rootOfSquareSum(phi[2].inner, phi0[2].inner);
				printers[0]->print("out.dat", mesh, phi);
				cout << step << "	" << error << endl;
			}
			step++;

		} while (step<maxStep&&error>convergenceThrehold);
		printers[0]->print("out.dat", mesh, phi);
	}


	void setRelativePressure(CellField &phi){
		for (int c = 0; c < mesh->cellNum; c++){
			phi.inner[c] = phi.inner[c] - phi.inner[0];
		}
	}
};









#endif