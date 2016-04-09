

#ifndef SOURCE_TERM_H_H
#define SOURCE_TERM_H_H

#include "serviceLocator.h"
#include "dataStructure.h"
class I_SourceTerm:public interface_t{

public :
	virtual vector<CSR> apply(vector<CellField>phi, Mesh* mesh, vector<double> s0) = 0;
	virtual vector<CSR> apply(vector<CellField*>phi, Mesh* mesh, vector<double> s0) = 0;
	virtual CSR apply(CellField phi, Mesh* mesh, vector<double> s0) = 0;

};


class SourceTerm :public I_SourceTerm, public member_t<SourceTerm>{

public:

	virtual vector<CSR>  apply(vector<CellField*>phi, Mesh* mesh, vector<double> s0) {
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return apply(phi1, mesh,s0);
	}

	virtual vector<CSR>  apply(vector<CellField>phi, Mesh* mesh, vector<double> s0) {
		vector<CSR> r;
		for (int i = 0; i < phi.size(); i++){
			vector<double> s0i(s0.begin() + i*mesh->cellNum, s0.begin() + (i + 1)*mesh->cellNum );
			r.push_back(apply(phi[i], mesh, s0i));
		}
		return r;
	}
	virtual CSR apply(CellField phi, Mesh* mesh,vector<double> s0){
		CSR eq(mesh);
		for (int c = 0; c<mesh->cellNum; c++){
			eq.b[c] -= mesh->C_v[c] * s0[c];
		}
		return eq;
	}

};






#endif