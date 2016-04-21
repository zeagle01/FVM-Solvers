#ifndef CONVECTION_H_H
#define CONVECTION_H_H


#include "Tool.h"
#include "faceReconstruct.h"







class UpwindConvection_Imp{
public:

	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return apply(V, phi1, mesh);
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
		vector<double> u_f(V.begin(), V.begin() + mesh->faceNum);
		vector<double> v_f(V.begin() + mesh->faceNum, V.begin() + mesh->faceNum * 2);
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

class DefferCenterConvection_Imp:public UpwindConvection_Imp{


	FaceReconstruct_BEB beb;
public :
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		CSR eq(mesh);
		vector<double> phi_f = beb.apply(phi, mesh);
		for (int f = 0; f < mesh->innerFaceNum; f++) {
			int inerF = mesh->IF[mesh->innerFaceBegin]+f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];

			double Sf_x = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 0];
			double Sf_y = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 1];
			double mdot = u_f[inerF] * Sf_x + v_f[inerF] * Sf_y;
			//cout << mdot << endl;
			double a_on = VectorMath<double>::myMin(mdot, 0);
			double a_oo = VectorMath<double>::myMax(mdot, 0);
			eq.A[mesh->ON[f * 2 + 0]] = a_on;
			eq.A[mesh->IA[owner]] += a_oo;

			//eq.b[owner] += a_oo*phi.inner[owner] + a_on*phi.inner[neighbor] + mdot*phi_f[inerF];//Deffer correction
			eq.b[owner] += a_oo*phi.inner[owner] + a_on*phi.inner[neighbor] - mdot*phi_f[inerF];//Deffer correction


			double a_no = VectorMath<double>::myMin(-mdot, 0);
			double a_nn = VectorMath<double>::myMax(-mdot, 0);
			eq.A[mesh->ON[f * 2 + 1]] = a_no;
			eq.A[mesh->IA[neighbor]] += a_nn;

			//eq.b[neighbor] += a_nn*phi.inner[neighbor] + a_no*phi.inner[owner] - mdot*phi_f[inerF];//Deffer correction
			eq.b[neighbor] += a_nn*phi.inner[neighbor] + a_no*phi.inner[owner] + mdot*phi_f[inerF];//Deffer correction
			//copy(eq->A.begin(), eq->A.end(), ostream_iterator<T>(cout, " "));

		}
		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];

				double Sf_x = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 0];
				double Sf_y = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 1];
				double mdot = u_f[boundaryF] * Sf_x + v_f[boundaryF] * Sf_y;

				double a_oo = VectorMath<double>::myMax(mdot, 0);
				double a_on = VectorMath<double>::myMin(mdot, 0);


				double a = phi.boundaryCondition[pat * 3 + 0];
				double b1 = phi.boundaryCondition[pat * 3 + 1];
				double c = phi.boundaryCondition[pat * 3 + 2];
				double d = mesh->F_d[boundaryF];

				eq.A[mesh->IA[owner]] += a_oo + a / d / (a / d + b1)*a_on;
				eq.b[owner] -= c / (a / d + b1)*a_on;
				//eq.b[owner] += a_oo*phi.inner[owner] + a_on*phi.boundary[f] + mdot*phi_f[boundaryF];
				eq.b[owner] += a_oo*phi.inner[owner] + a_on*phi.boundary[f] - mdot*phi_f[boundaryF];
			}
		}
		return eq;
	}
};

class ExplicitCenter_Imp :public UpwindConvection_Imp{


	FaceReconstruct_BEB beb;
public:
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		CSR eq(mesh);
		vector<double> phi_f = beb.apply(phi, mesh);
		for (int f = 0; f < mesh->innerFaceNum; f++) {
			int inerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];

			double Sf_x = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 0];
			double Sf_y = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 1];
			double mdot = u_f[inerF] * Sf_x + v_f[inerF] * Sf_y;
			//cout << mdot << endl;

			eq.b[owner] -=   mdot*phi_f[inerF];


			eq.b[neighbor] += mdot*phi_f[inerF];
			//copy(eq->A.begin(), eq->A.end(), ostream_iterator<T>(cout, " "));

		}
		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];

				double Sf_x = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 0];
				double Sf_y = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 1];
				double mdot = u_f[boundaryF] * Sf_x + v_f[boundaryF] * Sf_y;

				eq.b[owner] -=  mdot*phi_f[boundaryF];
			}
		}
		return eq;
	}
};

class ImplicitCenter_Imp :public UpwindConvection_Imp{


	FaceReconstruct_BEB beb;
public:
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		CSR eq(mesh);
		vector<double> phi_f = beb.apply(phi, mesh);
		for (int f = 0; f < mesh->innerFaceNum; f++) {
			int inerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[inerF * 2 + 0];
			int neighbor = mesh->FC[inerF * 2 + 1];

			double Sf_x = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 0];
			double Sf_y = mesh->F_a[inerF] * mesh->F_n[inerF * 2 + 1];
			double mdot = u_f[inerF] * Sf_x + v_f[inerF] * Sf_y;
			double x0 = mesh->C_c[owner * 2 + 0];
			double y0 = mesh->C_c[owner * 2 + 1];
			double x1 = mesh->F_c[inerF * 2 + 0];
			double y1 = mesh->F_c[inerF * 2 + 1];
			double x2 = mesh->C_c[neighbor * 2 + 0];
			double y2 = mesh->C_c[neighbor * 2 + 1];
			double g1 = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			double g2 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			//cout << mdot << endl;
			eq.A[mesh->ON[f * 2 + 0]] = g2 / (g1 + g2)*mdot;
			eq.A[mesh->IA[owner]] += g1 / (g1 + g2)*mdot;

			eq.A[mesh->ON[f * 2 + 1]] = g1 / (g1 + g2)*(-mdot);
			eq.A[mesh->IA[neighbor]] += g2 / (g1 + g2)*(-mdot);
			//copy(eq->A.begin(), eq->A.end(), ostream_iterator<double>(cout, " "));

		}
		for (int pat = 0; pat < mesh->boundaryPatchNum; pat++) {
			for (int f = mesh->IF[pat]; f < mesh->IF[pat + 1]; f++) {
				int boundaryF = f;
				int owner = mesh->FC[boundaryF * 2 + 0];

				double Sf_x = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 0];
				double Sf_y = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 1];
				double mdot = u_f[boundaryF] * Sf_x + v_f[boundaryF] * Sf_y;

				double convection_on = mdot;


				double a = phi.boundaryCondition[pat * 3 + 0];
				double b1 = phi.boundaryCondition[pat * 3 + 1];
				double c = phi.boundaryCondition[pat * 3 + 2];
				double d = mesh->F_d[boundaryF];

				eq.A[mesh->IA[owner]] +=  a / d / (a / d + b1)*convection_on;
				eq.b[owner] -= c / (a / d + b1)*convection_on;
			}
		}

		return eq;
	}
};












class I_Convection:public interface_t {
public:
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh) = 0;
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh) = 0;
	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh) = 0;
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh) = 0;
};

class UpwindConvection :public I_Convection, public member_t<UpwindConvection>{
public:
	UpwindConvection_Imp imp;
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		return imp.apply(V, phi, mesh);
	}

	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh){
		return imp.apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh){
		return  imp.apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		return  imp.apply(u_f,v_f, phi, mesh);
	}
};



class DefferCenterConvection :public I_Convection, public member_t<DefferCenterConvection>{
public:
	UpwindConvection_Imp* imp = new DefferCenterConvection_Imp();
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}

	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh){
		return  imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		return  imp->apply(u_f, v_f, phi, mesh);
	}

};

class ExplicitCenter :public I_Convection, public member_t<ExplicitCenter>{
public:
	UpwindConvection_Imp* imp = new ExplicitCenter_Imp();
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}

	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh){
		return  imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		return  imp->apply(u_f, v_f, phi, mesh);
	}

};


class ImplicitCenter :public I_Convection, public member_t<ImplicitCenter>{
public:
	UpwindConvection_Imp* imp = new ImplicitCenter_Imp();
	virtual vector<CSR>apply(vector<double> V, vector<CellField*> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}

	virtual vector<CSR>apply(vector<double> V, vector<CellField> phi, Mesh* mesh){
		return imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> V, CellField phi, Mesh* mesh){
		return  imp->apply(V, phi, mesh);
	}
	virtual CSR apply(vector<double> u_f, vector<double> v_f, CellField phi, Mesh* mesh){
		return  imp->apply(u_f, v_f, phi, mesh);
	}

};

#endif