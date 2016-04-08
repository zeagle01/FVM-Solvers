#ifndef GRADIENT_H_H
#define GRADIENT_H_H

#include <vector>
#include "dataStructure.h"
class I_GaussCellGradient :public interface_t {
public:
	virtual vector<CellField> apply(vector<double> faceValue, Mesh* mesh) = 0;

};


class GaussCellGradient :public I_GaussCellGradient, public member_t<GaussCellGradient>{
	virtual vector<CellField> apply(vector<double> faceValue, Mesh* mesh){

		vector<CellField> grad(2, CellField());
		grad[0].setUniformValue(mesh, 0);
		grad[1].setUniformValue(mesh, 0);

		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int ownner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			grad[0].inner[ownner] += mesh->F_n[innerF * 2 + 0] * mesh->F_a[innerF] / mesh->C_v[ownner] * faceValue[innerF];
			grad[1].inner[ownner] += mesh->F_n[innerF * 2 + 1] * mesh->F_a[innerF] / mesh->C_v[ownner] * faceValue[innerF];

			grad[0].inner[neighbor] -= mesh->F_n[innerF * 2 + 0] * mesh->F_a[innerF] / mesh->C_v[neighbor] * faceValue[innerF];
			grad[1].inner[neighbor] -= mesh->F_n[innerF * 2 + 1] * mesh->F_a[innerF] / mesh->C_v[neighbor] * faceValue[innerF];
		}

		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			int ownner = mesh->FC[boundaryF * 2 + 0];
			grad[0].inner[ownner] += mesh->F_n[boundaryF * 2 + 0] * mesh->F_a[boundaryF] / mesh->C_v[ownner] * faceValue[boundaryF];
			grad[1].inner[ownner] += mesh->F_n[boundaryF * 2 + 1] * mesh->F_a[boundaryF] / mesh->C_v[ownner] * faceValue[boundaryF];
		}
		return grad;
	}

};



class I_FaceGradient {
public:
	virtual vector<double> apply(CellField phi, Mesh* mesh) = 0;
};

class PlainFaceGradient :public  I_FaceGradient{

public:
	vector<double> apply(CellField phi, Mesh* mesh) {
		vector<double> grad(mesh->faceNum, 0);
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