#ifndef TEMPORAL_H_H
#define TEMPORAL_H_H



#include "dataStructure.h"

class I_Temporal:public interface_t{
public:
	double dt;
	virtual vector<CSR> apply(vector<CellField> phi, Mesh* mesh) = 0;
	virtual vector<CSR> apply(vector<CellField*> phi, Mesh* mesh) = 0;
	virtual CSR apply(CellField phi, Mesh* mesh) = 0;
};





class ImplictEuler :public I_Temporal, public member_t<ImplictEuler>{
public:

	virtual CSR apply(CellField phi, Mesh* mesh) {
		CSR eq(mesh);
		for (int c = 0; c<mesh->cellNum; c++){
			eq.A[mesh->DA[c]] = mesh->C_v[c] / dt;
			eq.b[c] = phi.inner[c] * mesh->C_v[c] / dt;
		}
		return eq;
	}
	virtual vector<CSR> apply(vector<CellField> phi, Mesh* mesh) {
		vector<CSR> r;
		for (int i = 0; i < phi.size(); i++){
			r.push_back(apply(phi[i], mesh));
		}
		return r;
	}

	virtual vector<CSR>  apply(vector<CellField*> phi, Mesh* mesh) {
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return apply(phi1, mesh);
	}
};


#endif