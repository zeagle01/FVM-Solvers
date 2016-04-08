#ifndef FACERECONSTRUCT_H_H
#define FACERECONSTRUCT_H_H

#include "dataStructure.h"
class I_FaceReconstruct :public interface_t{
public:
	virtual vector<double> apply(CellField phi, Mesh* mesh) = 0;
};

class FaceReconstruct_BEB :public I_FaceReconstruct, public member_t<FaceReconstruct_BEB>{
public:
	virtual vector<double> apply(CellField gama, Mesh* mesh) {
		vector<double> faceValue(mesh->faceNum);
		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int ownner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			double x0 = mesh->C_c[ownner * 2 + 0];
			double y0 = mesh->C_c[ownner * 2 + 1];
			double x1 = mesh->F_c[innerF * 2 + 0];
			double y1 = mesh->F_c[innerF * 2 + 1];
			double x2 = mesh->C_c[neighbor * 2 + 0];
			double y2 = mesh->C_c[neighbor * 2 + 1];
			double g1 = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			double g2 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			(faceValue)[innerF] = g1 / (g1 + g2)*gama.inner[neighbor] + g2 / (g1 + g2)*gama.inner[ownner];
		}
		boundaryHandle(gama, mesh, faceValue);
		return faceValue;
	}

	virtual void boundaryHandle(CellField gama, Mesh* mesh, vector<double> &faceValue) {
		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			faceValue[boundaryF] = gama.boundary[f];
		}
	}
};

class LinearFaceReconstruct_BEO :public I_FaceReconstruct, public member_t<LinearFaceReconstruct_BEO>{

public:
	virtual vector<double> apply(CellField gama, Mesh* mesh) {
		vector<double> faceValue(mesh->faceNum);
		for (int f = 0; f<mesh->innerFaceNum; f++){
			int innerF = mesh->IF[mesh->innerFaceBegin] + f;
			int ownner = mesh->FC[innerF * 2 + 0];
			int neighbor = mesh->FC[innerF * 2 + 1];
			double x0 = mesh->C_c[ownner * 2 + 0];
			double y0 = mesh->C_c[ownner * 2 + 1];
			double x1 = mesh->F_c[innerF * 2 + 0];
			double y1 = mesh->F_c[innerF * 2 + 1];
			double x2 = mesh->C_c[neighbor * 2 + 0];
			double y2 = mesh->C_c[neighbor * 2 + 1];
			double g1 = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			double g2 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
			(faceValue)[innerF] = g1 / (g1 + g2)*gama.inner[neighbor] + g2 / (g1 + g2)*gama.inner[ownner];
		}
		boundaryHandle(gama, mesh, faceValue);
		return faceValue;
	}
	virtual  void boundaryHandle(CellField gama, Mesh* mesh, vector<double> &faceValue) {
		for (int f = 0; f<mesh->boundaryFaceNum; f++){
			int boundaryF = f;
			int owner = mesh->FC[boundaryF * 2 + 0];
			faceValue[boundaryF] = gama.inner[owner];
		}
	}
};


/*
class B_RhieChow{
LinearFaceReconstruct_BEO beo;
FaceReconstruct_BEB beb;
public:
vector<CellField> calculateD(CSR u_momentum, CSR v_momentum, Mesh* mesh){
vector<CellField> D_f(mesh->cellNum * 2);
vector<CellField> D(2, CellField());
D[0].setUniformValue(mesh, 0);
D[1].setUniformValue(mesh, 0);
for (int c = 0; c<mesh->cellNum; c++){
D[0].inner[c] = mesh->C_v[c] / u_momentum.A[mesh->IA[c]];
D[1].inner[c] = mesh->C_v[c] / v_momentum.A[mesh->IA[c]];
}

vector<double> df0=beo.apply(D[0], mesh);
D_f.assign(df0.begin(), df0.end());
D_f[1] = beo.apply(D[1], mesh);
return D_f;
}
vector<vector<double>> reconstructVelocity(vector<vector<double>> D_f, vector< CellField>& p_cellGrad, vector<vector<double>> p_faceGrad, vector<CellField> velocity, Mesh* mesh){

vector<vector<double>> v_RC(2, vector<double>());
v_RC[0].resize(mesh->faceNum);
v_RC[1].resize(mesh->faceNum);
vector<double> p_avrgFaceGrad_x = beo.apply(p_cellGrad[0], mesh);
vector<double> p_avrgFaceGrad_y = beo.apply(p_cellGrad[1], mesh);
vector<double> u_avrg = beb.apply(velocity[0], mesh);
vector<double> v_avrg = beb.apply(velocity[1], mesh);
for (int f = 0; f<mesh->faceNum; f++){
v_RC[0][f] = u_avrg[f] + D_f[0][f] * (p_avrgFaceGrad_x[f] - p_faceGrad[0][f]);
v_RC[1][f] = v_avrg[f] + D_f[1][f] * (p_avrgFaceGrad_y[f] - p_faceGrad[1][f]);
}
return v_RC;
}


};
*/

#endif