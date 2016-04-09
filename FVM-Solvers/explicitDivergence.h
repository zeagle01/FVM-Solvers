

#ifndef EXPLICIT_DIVERGENCE_H_H
#define EXPLICIT_DIVERGENCE_H_H


#include "dataStructure.h"
class I_ExplicitDivergence :public interface_t{

public:
	virtual vector<double> apply(vector<double> V_f, Mesh* mesh) = 0;

};


class ExplicitDivergence :public I_ExplicitDivergence, public member_t<ExplicitDivergence>{

public:
	virtual vector<double> apply(vector<double> V_f, Mesh* mesh){
		vector<double> div(mesh->cellNum);
		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int ownner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			double face_x = mesh->F_a[innerF] * mesh->F_n[innerF * 2 + 0];
			double face_y = mesh->F_a[innerF] * mesh->F_n[innerF * 2 + 1];

			double flux = face_x*V_f[innerF] + face_y*V_f[mesh->faceNum + innerF];
			div[ownner] += flux / mesh->C_v[ownner];
			div[neighbor] -= flux / mesh->C_v[neighbor];
		}
		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			int ownner = mesh->FC[boundaryF * 2 + 0];
			double face_x = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 0];
			double face_y = mesh->F_a[boundaryF] * mesh->F_n[boundaryF * 2 + 1];
			double flux = face_x*V_f[boundaryF] + face_y*V_f[boundaryF + mesh->faceNum];
			div[ownner] += flux / mesh->C_v[ownner];
		}
		return div;
	}
};

#endif