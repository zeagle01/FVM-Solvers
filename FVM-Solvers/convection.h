#ifndef CONVECTION_H_H
#define CONVECTION_H_H


#include "Tool.h"
#include "faceReconstruct.h"
class I_Convection:public interface_t {
public:
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh) = 0;
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh) = 0;
	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh) = 0;
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh) = 0;
};

class UpwindConvection :public I_Convection, public member_t<UpwindConvection>{
public:

	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return apply(V, phi1,mesh);
	}

	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh){
		vector<CSR> r;
		for (int i = 0; i < phi.size(); i++){
			CSR eq = apply(V, phi[0], mesh);
			r.push_back(eq);
		}
		return r;
	}
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh){
		CSR eq(mesh);
		vector<double> u_f(V.begin(), V.begin() + mesh->faceNum );
		vector<double> v_f( V.begin() + mesh->faceNum,V.begin()+mesh->faceNum*2);
		return apply(u_f, v_f, phi, mesh);
	}
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		CSR eq(mesh);
		//copy(eq->A.begin(), eq->A.end(), ostream_iterator<double>(cout, " "));

		for (int f = 0; f < mesh->innerFaceNum; f++) {
			int inerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];

			double Sf_x = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 0];
			double Sf_y = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 1];
			double mdot = u_f[inerF] * Sf_x + v_f[inerF] * Sf_y;
			//cout << mdot << endl;
			eq.A[mesh->ON[f * 2 + 0]] = VectorMath<double>::myMin(mdot, 0);
			eq.A[mesh->IA[owner]] += VectorMath<double>::myMax(mdot, 0);

			eq.A[mesh->ON[f * 2 + 1]] = VectorMath<double>::myMin(-mdot, 0);
			eq.A[mesh->IA[neighbor]] += VectorMath<double>::myMax(-mdot, 0);
			//copy(eq->A.begin(), eq->A.end(), ostream_iterator<double>(cout, " "));

		}
		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];

				double Sf_x = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 0];
				double Sf_y = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 1];
				double mdot = u_f[boundaryF] * Sf_x + v_f[boundaryF] * Sf_y;

				double convection_oo = VectorMath<double>::myMax(mdot, 0);
				double convection_on = VectorMath<double>::myMin(mdot, 0);


				double a = phi.boundaryCondition[pat * 3 + 0];
				double b1 = phi.boundaryCondition[pat * 3 + 1];
				double c = phi.boundaryCondition[pat * 3 + 2];
				double d = mesh->F_d[boundaryF];

				eq.A[mesh->IA[owner]] += convection_oo + a / d / (a / d + b1)*convection_on;
				eq.b[owner] -= c / (a / d + b1)*convection_on;
			}
		}

		return eq;
	}
};


#endif