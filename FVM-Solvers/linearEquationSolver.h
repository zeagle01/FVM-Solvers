#ifndef LINEAREQUATIONSOLVER_H_H
#define LINEAREQUATIONSOLVER_H_H


#include "dataStructure.h"
#include "mesh.h"
#include "serviceLocator.h"
class I_LinearEquationSolver :public interface_t{

public:
	double converge_threhold;
	int max_step;
	int check_step;
	virtual void solve(CSR eq, CellField& phi, Mesh* mesh){
		double error = 0; int step = 0;
		CellField pre;
		do{
			if (step%check_step == 0){
				pre = phi;
			}
			iterate(eq, phi, mesh, 1);
			if (step%check_step == 0){
				error = VectorMath<double>::rootOfSquareSum(phi.inner.data(), pre.inner.data(), phi.inner.size());
			}
			step++;
		} while (error>converge_threhold&&step<max_step);
	}
	virtual void iterate(CSR eq, CellField& phi, Mesh* mesh, int iterateTimes) = 0;
};


class  GS_Solver : public I_LinearEquationSolver, public member_t<GS_Solver>{
public:
	virtual void iterate(CSR eq, CellField& phi, Mesh* mesh, int iterateTimes) {
		for (int it = 0; it < iterateTimes; ++it){
			for (int c = 0; c < mesh->cellNum; c++) {
				double temp = 0;
				for (int i = mesh->IA[c]; i < mesh->IA[c + 1]; i++) {
					if (mesh->JA[i] != c)
					temp += eq.A[i] * phi.inner[mesh->JA[i]];
				}
				phi.inner[c] = (eq.b[c] - temp) / eq.A[mesh->DA[c]];
			}
		}
	}

};



class GMRES : public I_LinearEquationSolver, public member_t<GMRES>{

public:
	int m = 5;
	virtual void iterate(CSR eq, CellField& phi, Mesh* mesh, int iterateTimes) {

		
		if (VectorMath<double>::length(eq.b.data(), mesh->cellNum)<1e-50&&VectorMath<double>::length(phi.inner.data(), mesh->cellNum)<1e-50){
			phi.inner[0] += 1.0;
		}
		vector<double> x = phi.inner;
		vector<double> r = x;
		for (int it = 0; it<iterateTimes; it++){
			CSR::MatrixVectorMul(eq.A.data(), mesh->IA.data(), mesh->JA.data(),x.data(),mesh->cellNum, r.data());
			VectorMath<double>::plus_minus(-1, eq.b.data(), mesh->cellNum,r.data());
			VectorMath<double>::reverse(r.data(), mesh->cellNum);
			double beta = VectorMath<double>::length(r.data(), mesh->cellNum);
			vector<vector<double>> V (m+1,vector<double>());
			
			VectorMath<double>::scalarMul(1.0 / beta, mesh->cellNum, r.data());
			V[0] = r;
			vector<vector<double>> H (m+1,vector<double>());
			for (int i = 0; i < m + 1; i++){
				H[i].resize(m);
			}
			vector<vector<double>> R(m, vector<double>());
			for (int i = 0; i < m ; i++){
				R[i].resize(m);
			}
			vector<double> c(m);
			vector<double> s(m);
			vector<double> b_hat(m + 1); 
			vector<double> y(m);
			b_hat[0] = beta;
			int i = 0;
			for (i = 0; i<m; i++){
				vector<double> w(mesh->cellNum);
				CSR::MatrixVectorMul(eq.A.data(), mesh->IA.data(), mesh->JA.data(), V[i].data(),mesh->cellNum,w.data());
				for (int j = 0; j <= i; j++){
					H[j][i] = VectorMath<double>::dot(w.data(), V[j].data(),mesh->cellNum);
					VectorMath<double>::plus_minus( -H[j][i], V[j].data(),mesh->cellNum,w.data());
				}
				H[i + 1][i] = VectorMath<double>::length(w.data(), mesh->cellNum);
				VectorMath<double>::scalarMul(1.0 / H[i + 1][i],mesh->cellNum, w.data());
				V[i + 1] = w;
				R[0][i] = H[0][i];
				for (int j = 1; j <= i; j++){
					double gama = c[j - 1] * R[j - 1][i] + s[j - 1] * H[j][i];
					R[j][i] = -s[j - 1] * R[j - 1][i] + c[j - 1] * H[j][i];
					R[j - 1][i] = gama;
				}
				double delta =sqrt(R[i][i] * R[i][i] + H[i + 1][i] * H[i + 1][i]);
				c[i] = R[i][i] / delta;
				s[i] = H[i + 1][i] / delta;
				R[i][i] = c[i] * R[i][i] + s[i] * H[i + 1][i];
				b_hat[i + 1] = -s[i] * b_hat[i];
				b_hat[i] = c[i] * b_hat[i];
				double rho = abs(b_hat[i + 1]);
				if (rho<converge_threhold){
					i++;
					break;
				}
			}

			y[i - 1] = b_hat[i - 1] / R[i - 1][i - 1];
			for (int k = i - 2; k >= 0; k--){
				double temp = 0;
				for (int l = k + 1; l<i; l++){
					temp += R[k][l] * y[l];
				}
				y[k] = (b_hat[k] - temp) / R[k][k];
			}
			for (int k = 0; k<mesh->cellNum; k++){
				for (int l = 0; l<i; l++){
					x[k] += y[l] * V[l][k];
				}
			}
		}
		phi.inner = x;
	}
};

/*
class  OpenMP_Jacobi_Solver : public I_LinearEquationSolver, public member_t< OpenMP_Jacobi_Solver>{
public:
virtual void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes) {
for (int it = 0; it < iterateTimes; ++it){
vector<double> temp(mesh.cellNum, 0);
#pragma omp parallel for
for (int c = 0; c < mesh.cellNum; c++) {
for (int i = mesh.IA[c] + 1; i < mesh.IA[c + 1]; i++) {
temp[c] += eq.A[i] * phi.inner[mesh.JA[i]];
}
phi.inner[c] = (eq.b[c] - temp[c]) / eq.A[mesh.IA[c]];
}
}
}

};


class  GPU_Jacobi_Solver : public I_LinearEquationSolver, public member_t<GPU_Jacobi_Solver>{
public:
virtual void solve(CSR eq, CellField& phi, Mesh mesh);
void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes){};
};
*/

/*
class  GPU_Jacobi_Solver_vector : public I_LinearEquationSolver, public member_t<GPU_Jacobi_Solver_vector>{


public:
virtual void solve(CSR eq, CellField& phi, Mesh mesh);
void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes);
};
*/

#endif