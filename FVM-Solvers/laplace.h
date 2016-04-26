#ifndef LAPLACE_H_H
#define LAPLACE_H_H

#include "dataStructure.h"
#include <vector>
#include "serviceLocator.h"
class I_Laplace :public interface_t{
public:
	vector<double> grad_f;
	virtual CSR apply(vector<double> gama_f, CellField phi, Mesh* mesh) = 0;
	virtual CSR apply(double  gama_f, CellField phi, Mesh* mesh) = 0;
	virtual vector<CSR> apply(double gama_f, vector<CellField> phi, Mesh* mesh) = 0;
	virtual vector<CSR> apply(double gama_f, vector<CellField*> phi, Mesh* mesh) = 0;
};



class Laplace :public I_Laplace, public member_t<Laplace> {
public:

	virtual vector<CSR> apply(double gama_f, vector<CellField*> phi, Mesh* mesh){
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return apply(gama_f, phi1, mesh);
	}

	virtual vector<CSR> apply(double gama_f, vector<CellField> phi, Mesh* mesh){
		vector<CSR> r;
		vector<double> gama_field(mesh->faceNum, gama_f);
		for (int i = 0; i < phi.size(); i++){

			CSR eq = apply(gama_field, phi[i], mesh);
			r.push_back(eq);
		}
		return r;
	}


	virtual CSR apply(double gama_f, CellField phi, Mesh* mesh){
		vector<double> gama_field(mesh->faceNum, gama_f);
		return apply(gama_field, phi, mesh);
	}

	virtual CSR apply(vector<double> gama_f, CellField phi, Mesh* mesh){
		CSR eq(mesh);


		for (int f = 0; f < mesh->innerFaceNum; f++){
			int inerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];
			double S_f = mesh->F_a[inerF];
			double d = mesh->F_d[inerF];
			double diffusion_on = gama_f[inerF] * S_f / d;
			double diffusion_oo = -diffusion_on;
			double diffusion_no = diffusion_on;
			double diffusion_nn = -diffusion_on;

			eq.A[mesh->DA[owner]] += diffusion_oo;
			eq.A[mesh->ON[f * 2 + 0]] += diffusion_on;
			eq.A[mesh->DA[neighbor]] += diffusion_nn;
			eq.A[mesh->ON[f * 2+1]] += diffusion_no;
		}

		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];
				double S_f = mesh->F_a[boundaryF];
				double d = mesh->F_d[boundaryF];

				double diffusion_on = gama_f[boundaryF] * S_f / d;
				double diffusion_oo = -diffusion_on;
				double a = phi.boundaryCondition[pat * 3 + 0];
				double b1 = phi.boundaryCondition[pat * 3 + 1];
				double c = phi.boundaryCondition[pat * 3 + 2];
				eq.A[mesh->DA[owner]] += diffusion_oo + a / d / (a / d + b1)*diffusion_on;
				eq.b[owner] -= c / (a / d + b1)*diffusion_on;
			}
		}
		return eq;
	}

};



class UnorthogonalCorrectionLaplace :public I_Laplace, public member_t<UnorthogonalCorrectionLaplace> {
public:
	Laplace laplace;
	virtual CSR apply(double gama_f, CellField phi, Mesh* mesh){
		vector<double> gama_field(mesh->faceNum, gama_f);
		return apply(gama_field, phi, mesh);
	}

	virtual CSR apply(vector<double> gama_f, CellField phi, Mesh* mesh){
		CSR eq = laplace.apply(gama_f, phi, mesh);

		for (int f = 0; f < mesh->innerFaceNum; f++){
			int inerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];

			double t0 = mesh->F_n[inerF * 2 + 0] - mesh->F_e[inerF * 2 + 0];
			double t1 = mesh->F_n[inerF * 2 + 1] - mesh->F_e[inerF * 2 + 1];

			double S_f = mesh->F_a[inerF];

			double cf = S_f*(t0*grad_f[inerF] + t1*grad_f[inerF + mesh->faceNum]);
			eq.b[owner] -= cf;
			eq.b[neighbor] += cf;
		}

		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];
				double t0 = mesh->F_n[boundaryF * 2 + 0] - mesh->F_e[boundaryF * 2 + 0];
				double t1 = mesh->F_n[boundaryF * 2 + 1] - mesh->F_e[boundaryF * 2 + 1];

				double S_f = mesh->F_a[boundaryF];

				double cf = S_f*(t0*grad_f[boundaryF] + t1*grad_f[boundaryF + mesh->faceNum]);
				eq.b[owner] -= cf;

				double a = phi.boundaryCondition[pat * 3 + 0];
				eq.b[owner] -= a*cf;
			}
		}
		return eq;
	}

};



#endif