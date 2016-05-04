#ifndef FACERECONSTRUCT_H_H
#define FACERECONSTRUCT_H_H

#include "dataStructure.h"
class I_FaceReconstruct :public interface_t{
public:
	virtual vector<double> apply(CellField phi, Mesh* mesh) = 0;
	virtual vector<double> apply(vector<CellField> phi, Mesh* mesh) = 0;
	virtual vector<double> apply(vector<CellField*> phi, Mesh* mesh) = 0;
};

class FaceReconstruct_BEB :public I_FaceReconstruct, public member_t<FaceReconstruct_BEB>{
public:

	virtual vector<double> apply(vector<CellField> phi, Mesh* mesh){
		vector<double> r;
		for (int i = 0; i < phi.size(); i++){
			vector<double> fv = apply(phi[i], mesh);
			r.insert(r.end(), fv. begin(), fv.end());
		}
		return r;
	}
	virtual vector<double> apply(vector<CellField*> phi, Mesh* mesh){
		vector<double> r;
		for (int i = 0; i < phi.size(); i++){
			vector<double> fv = apply(*(phi[i]), mesh);
			r.insert(r.end(), fv.begin(), fv.end());
		}
		return r;
	}

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

	virtual vector<double> apply(vector<CellField> phi, Mesh* mesh){
		vector<double> r;
		for (int i = 0; i < phi.size(); i++){
			vector<double> fv = apply(phi[i], mesh);
			r.insert(r.end(), fv.begin(), fv.end());
		}
		return r;
	}
	virtual vector<double> apply(vector<CellField*> phi, Mesh* mesh){
		vector<double> r;
		for (int i = 0; i < phi.size(); i++){
			vector<double> fv = apply(*(phi[i]), mesh);
			r.insert(r.end(), fv.begin(), fv.end());
		}
		return r;
	}

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



class B_RhieChow:public interface_t{

public:
	vector<I_FaceReconstruct*> faceReconstruct;
	virtual vector<double> calculateD(vector<CSR> V_eq, Mesh* mesh) = 0;
		
	virtual vector<double> reconstructVelocity(vector<double> D_f, vector<CellField> p_cellGrad, vector<double> p_faceGrad, vector<CellField> velocity, Mesh* mesh) = 0;
	virtual vector<double> reconstructVelocity(vector<double> D_f, vector<CellField> p_cellGrad, vector<double> p_faceGrad, vector<CellField*> velocity, Mesh* mesh) = 0;
		

};




class RhieChowInterpolation :public B_RhieChow, public member_t<RhieChowInterpolation> {

public:
	virtual vector<double> calculateD(vector<CSR> V_eq, Mesh* mesh) {
		vector<double> D_f(mesh->cellNum * 2);
		vector<CellField> D(2, CellField());
		D[0].setUniformValue(mesh, 0);
		D[1].setUniformValue(mesh, 0);
		for (int c = 0; c < mesh->cellNum; c++){
			D[0].inner[c] = mesh->C_v[c] / V_eq[0].A[mesh->DA[c]];
			D[1].inner[c] = mesh->C_v[c] / V_eq[1].A[mesh->DA[c]];
		}
		vector<double> df0 = faceReconstruct[1]->apply(D[0], mesh);
		vector<double> df1 = faceReconstruct[1]->apply(D[1], mesh);
		D_f.assign(df0.begin(), df0.end());
		D_f.insert(D_f.end(), df1.begin(), df1.end());
		return D_f;
	}

	virtual vector<double> reconstructVelocity(vector<double> D_f, vector<CellField> p_cellGrad, vector<double> p_faceGrad, vector<CellField> velocity, Mesh* mesh){
		vector<double> v_RC(mesh->faceNum*2);
		vector<double> p_avrgFaceGrad = faceReconstruct[1]->apply(p_cellGrad, mesh);
		vector<double> V_avrg = faceReconstruct[0]->apply(velocity, mesh);
		for (int f = 0; f < mesh->faceNum; f++){
			v_RC[f] = V_avrg[f] + D_f[f] * (p_avrgFaceGrad[f] - p_faceGrad[f]);
			v_RC[mesh->faceNum + f] = V_avrg[mesh->faceNum + f] + D_f[mesh->faceNum + f] * (p_avrgFaceGrad[mesh->faceNum + f] - p_faceGrad[mesh->faceNum + f]);
		}
		return v_RC;
	}
	virtual vector<double> reconstructVelocity(vector<double> D_f, vector<CellField> p_cellGrad, vector<double> p_faceGrad, vector<CellField*> phi, Mesh* mesh){
		vector<CellField> phi1;
		for (int i = 0; i < phi.size(); i++){
			phi1.push_back(*(phi[i]));
		}
		return reconstructVelocity( D_f,  p_cellGrad,p_faceGrad,  phi1,  mesh);
	}
};


#endif