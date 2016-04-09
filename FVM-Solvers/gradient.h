#ifndef GRADIENT_H_H
#define GRADIENT_H_H

#include <vector>
#include "dataStructure.h"
class I_GaussCellGradient :public interface_t {
public:
	virtual vector<double> apply(vector<double> faceValue, Mesh* mesh) = 0;

};


class GaussCellGradient :public I_GaussCellGradient, public member_t<GaussCellGradient>{
	virtual vector<double> apply(vector<double> faceValue, Mesh* mesh){

		vector<double> grad(2*mesh->cellNum);

		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int ownner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			grad[ownner] += mesh->F_n[innerF * 2 + 0] * mesh->F_a[innerF] / mesh->C_v[ownner] * faceValue[innerF];
			grad[ownner+mesh->cellNum] += mesh->F_n[innerF * 2 + 1] * mesh->F_a[innerF] / mesh->C_v[ownner] * faceValue[innerF];

			grad[neighbor] -= mesh->F_n[innerF * 2 + 0] * mesh->F_a[innerF] / mesh->C_v[neighbor] * faceValue[innerF];
			grad[neighbor+mesh->cellNum] -= mesh->F_n[innerF * 2 + 1] * mesh->F_a[innerF] / mesh->C_v[neighbor] * faceValue[innerF];
		}

		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			int ownner = mesh->FC[boundaryF * 2 + 0];
			grad[ownner] += mesh->F_n[boundaryF * 2 + 0] * mesh->F_a[boundaryF] / mesh->C_v[ownner] * faceValue[boundaryF];
			grad[ownner+mesh->cellNum] += mesh->F_n[boundaryF * 2 + 1] * mesh->F_a[boundaryF] / mesh->C_v[ownner] * faceValue[boundaryF];
		}
		return grad;
	}

};



class I_FaceGradient :public interface_t{
public:
	virtual vector<double> apply(CellField phi, Mesh* mesh) = 0;
};

class PlainFaceGradient :public  I_FaceGradient, public member_t<PlainFaceGradient>{

public:
	vector<double> apply(CellField phi, Mesh* mesh) {
		vector<double> grad(mesh->faceNum*2, 0);
		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int owner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			double phi2 = phi.inner[neighbor];
			double phi1 = phi.inner[owner];
			double temp = (phi2 - phi1) / mesh->F_d[innerF];
			grad[innerF] = temp*mesh->F_n[innerF * 2 + 0];
			grad[mesh->faceNum + innerF] = temp*mesh->F_n[innerF * 2 + 1];
		}
		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			int owner = mesh->FC[boundaryF * 2 + 0];
			double phi2 = phi.boundary[f];
			double phi1 = phi.inner[owner];
			double temp = (phi2 - phi1) / mesh->F_d[boundaryF];
			grad[boundaryF] = temp*mesh->F_n[boundaryF * 2 + 0];
			grad[mesh->faceNum + boundaryF] = temp*mesh->F_n[boundaryF * 2 + 1];
		}
		return grad;
	}

};


#endif