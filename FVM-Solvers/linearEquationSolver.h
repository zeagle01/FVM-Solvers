#ifndef LINEAREQUATIONSOLVER_H_H
#define LINEAREQUATIONSOLVER_H_H


#include "dataStructure.h"
#include "mesh.h"
#include "serviceLocator.h"




class I_Preconditioner :public interface_t{
public:
	virtual void solve(double* A, int*DA, int* IA, int* JA, double* b, int n, double* x) = 0;
};

class GS_Precondition : public I_Preconditioner, public member_t<GS_Precondition>{
public:
	virtual void solve(double* A, int*DA, int* IA, int* JA, double* b, int n, double* x) {
		x[0] = b[0] / A[0];
		for (int i = 1; i < n; i++){
			double temp = 0;
			int j = 0;
			for (j = IA[i]; j < IA[i + 1]; j++){
				if (JA[j] < i){
					temp += A[j] * x[JA[j]];
				}
				else if (JA[j] == i){
					break;
				}
			}
			x[i] = (b[i] - temp) / A[j];
		}
	}
};

class Jacobi_Precondition : public I_Preconditioner, public member_t<Jacobi_Precondition>{
public:
	virtual void solve(double* A, int*DA, int* IA, int* JA, double* b, int n, double* x) {
		for (int i = 0; i < n; i++){
			for (int j = IA[i]; j < IA[i + 1]; j++){
				if (JA[j] == i){
					x[i] = b[i] / A[j];
				}
			}
		}
	}
};

class LUSGS_Precondition : public I_Preconditioner, public member_t<LUSGS_Precondition>{
public:
	virtual void solve(double* A, int*DA, int* IA, int* JA, double* b, int n, double* x) {
		vector<double> y_v(n);
		double *y = y_v.data();
		y[0] = b[0];
		for (int i = 1; i < n; i++){
			double temp = 0;
			int j = 0;
			for (j = IA[i]; j < IA[i + 1]; j++){
				if (JA[j] < i){
					temp += A[j] / A[DA[JA[j]]] * y[JA[j]];
				}
				else if (JA[j] == i){
					break;
				}
			}
			y[i] = (b[i] - temp);
		}
		x[n - 1] = y[n - 1] / A[DA[n - 1]];
		for (int i = n - 2; i >= 0; i--){
			double temp = 0;
			for (int j = IA[i]; j < IA[i + 1]; j++){
				if (JA[j]>i){
					temp += A[j] * x[JA[j]];
				}
				else if (JA[j] == i){
					break;
				}
			}
			x[i] = (y[i] - temp) / A[DA[i]];
		}
	}
};


class I_LinearEquationSolver :public interface_t{

public:
	double converge_threhold;
	int max_step;
	int check_step;
	I_Preconditioner* pre;
	virtual void solve(CSR eq, CellField& phi, Mesh* mesh) = 0;
};


class  GS_Solver : public I_LinearEquationSolver, public member_t<GS_Solver>{
public:
	virtual void solve(CSR eq, CellField& phi, Mesh* mesh){
		for (int it = 0; it < max_step; ++it){
			vector<double> x_old = phi.inner;
			for (int c = 0; c < mesh->cellNum; c++) {
				double temp = 0;
				for (int i = mesh->IA[c]; i < mesh->IA[c + 1]; i++) {
					if (mesh->JA[i] != c)
						temp += eq.A[i] * phi.inner[mesh->JA[i]];
				}
				phi.inner[c] = (eq.b[c] - temp) / eq.A[mesh->DA[c]];
			}
			
			if (it%check_step == 0){
				double e_x = VectorMath<double>::rootOfSquareSum(phi.inner.data(), x_old.data(), mesh->cellNum);
				if (e_x < converge_threhold) {	
					break;
				}
			}
		}
	}

};

class  Jacobi_Solver : public I_LinearEquationSolver, public member_t<Jacobi_Solver>{
public:
	virtual void solve(CSR eq, CellField& phi, Mesh* mesh){
		for (int it = 0; it < max_step; ++it){
			vector<double> x_old = phi.inner;
			for (int c = 0; c < mesh->cellNum; c++) {
				double temp = 0;
				for (int i = mesh->IA[c]; i < mesh->IA[c + 1]; i++) {
					if (mesh->JA[i] != c)
						temp += eq.A[i] * x_old[mesh->JA[i]];
				}
				phi.inner[c] = (eq.b[c] - temp) / eq.A[mesh->DA[c]];
			}

			if (it%check_step == 0){
				double e_x = VectorMath<double>::rootOfSquareSum(phi.inner.data(), x_old.data(), mesh->cellNum);
				if (e_x < converge_threhold) {
					break;
				}
			}
		}
	}
	//virtual void iterate(CSR eq, CellField& phi, Mesh* mesh, int iterateTimes) {}

};








class GMRES : public I_LinearEquationSolver, public member_t<GMRES>{

public:
	int m = 5;
	virtual void solve(CSR eq, CellField& phi, Mesh* mesh) {
		double* x = phi.inner.data();
		int n = mesh->cellNum;
		double* A = eq.A.data();
		int* IA = mesh->IA.data();
		int* JA = mesh->JA.data();
		double* b = eq.b.data();


		for (int it = 0; it < max_step; it++){
			vector <double> x_old_v(n);
			double * x_old = x_old_v.data();
			vector<double> V_v(n*(m + 1));
			double *V = V_v.data();
			vector<double> H_v(m*(m + 1));
			double* H = H_v.data();
			vector<double> R_v(m*m);
			double* R = R_v.data();
			vector<double> c_v(m);
			vector<double> s_v(m);
			double* c = c_v.data();
			double* s = s_v.data();
			vector<double> b_hat_v(m + 1);
			double* b_hat = b_hat_v.data();
			vector<double> y_v(m);
			double* y = y_v.data();

			memcpy(x_old, x, n*sizeof(double));
			CSR::MatrixVectorMul(A, IA, JA, x, n,V);
			VectorMath<double>::plus_minus(-1, b, n, V);
			VectorMath<double>::reverse(V,n);
			if (pre){
				vector<double> temp_v(V, V + n);
				double* temp = temp_v.data();
				pre->solve(A, mesh->DA.data(), IA, JA, temp, n, V);
			}
			double beta = VectorMath<double>::length(V, n);
			b_hat[0] = beta;
			VectorMath<double>::scalarMul(1.0 / beta, n, V);
			int i = 0;
			for (i = 0; i < m; i++){
				CSR::MatrixVectorMul(A, IA, JA, V + i*n, n,V + (i + 1)*n);
				if (pre){
					vector<double> temp_v(V + (i + 1)*n, V + (i + 2)*n);
					double* temp = temp_v.data();
					pre->solve(A, mesh->DA.data(), IA, JA, temp, n, V + (i + 1)*n);
				}
				for (int j = 0; j <= i; j++){
					H[j + i*(m + 1)] = VectorMath<double>::dot(V + j*n, V + (i + 1)*n, n);
					VectorMath<double>::plus_minus(-H[j + i*(m + 1)], V + j*n, n, V + (i + 1)*n);
				}
				H[i + 1 + i*(m + 1)] = VectorMath<double>::length(V + (i + 1)*n, n);
				VectorMath<double>::scalarMul(1.0 / H[i + 1 + i*(m + 1)], n, V + (i + 1)*n);


				R[0 + i*m] = H[0 + i*(m + 1)];
				for (int j = 0; j < i; j++){
					double temp = c[j] * R[j + i*m] + s[j] * H[j + 1 + i*(m + 1)];
					R[j+1+i*m] = -s[j] * R[j + i*m] + c[j] * H[j + 1 + i*(m + 1)];
					R[j + i*m] = temp;
				}
				double rr = R[i + i*m];
				double hh = H[i + 1 + i*(m + 1)];
				double delta = sqrt(rr*rr + hh*hh);
				c[i] = rr / delta;
				s[i] = hh / delta;
				R[i + i*m] = c[i] * R[i + i*m] + s[i] * H[i + 1 + i*(m + 1)];
				b_hat[i + 1] = -s[i] * b_hat[i];
				b_hat[i] = c[i] * b_hat[i];
				double e_r = abs(b_hat[i + 1]);
				if (e_r < converge_threhold){
					i++;
					break;
				}
			}
			
			y[i - 1] = b_hat[i - 1] / R[i - 1 + (i - 1)*m];
			for (int k = i - 2; k >= 0; k--){
				double temp = 0;
				for (int l = k + 1; l < i; l++){
					temp += R[k + l*m] * y[l];
				}
				y[k] = (b_hat[k] - temp) / R[k + k*m];
			}
			
			for (int k = 0; k < i; k++){
				VectorMath<double>::plus_minus(y[k], V + k*n, n,x);
			}
			if (it%check_step == 0){
				double e_x = VectorMath<double>::rootOfSquareSum(x, x_old, n);
				if (e_x < converge_threhold) {
					break;
				}
			}
		}
	}
	//virtual void iterate(CSR eq, CellField& phi, Mesh* mesh, int iterateTimes) {}
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