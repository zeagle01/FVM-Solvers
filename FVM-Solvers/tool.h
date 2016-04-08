#ifndef TOOL_H_H
#define TOOL_H_H


#include <vector>

template<class T>class VectorMath {
public:
	static T outerdot(T x0, T y0, T x1, T y1) {
		return x0*y1 - x1*y0;
	}
	static T triangleArea(T dx0, T dy0, T dx1, T dy1) {
		return 0.5*abs(dx0*dy1 - dx1*dy0);
	}

	static T* coorthogonal(T dx, T dy) {
		T* v = new T[2];
		//T[2] v ;
		v[0] = dy;
		v[1] = -dx;
		return v;
	}
	static T dot(T x0, T y0, T x1, T y1) {

		return x0*x1 + y0*y1;
	}
	static T* normalize(T dx, T dy) {
		T* norm = new T[2];
		//T[2] norm;
		T l = length(dx, dy);
		norm[0] = dx / l;
		norm[1] = dy / l;
		return norm;
	}
	static T length(T dx, T dy) {
		return sqrt(dx*dx + dy*dy);
	}
	static T* reverse(T dx, T dy) {
		T* r = new T[2];
		//T[2] r ;
		r[0] = -dx;
		r[1] = -dy;
		return r;
	}

	static T myMax(T a, T b){
		return a > b ? a : b;
	}

	static T myMin(T a, T b){

		return a < b ? a : b;
	}

	/*
	static T rootOfSquareSum(vector<T> E1, vector<T> E2){
	transform(E1.cbegin(), E1.cend(), E2.cbegin(), E2.begin(), minus<T>());
	transform(E2.begin(), E2.end(), E2.begin(), E2.begin(), multiplies<T>());
	T error = accumulate(E2.begin(), E2.end(), (T)0);//have to force convert to T type unless the result is always 0
	//error = inner_product(pre.inner.cbegin(), pre.inner.cend(), pre.inner.cbegin(),
	//	1, plus<T>(), plus<T>());
	return sqrt(error);
	}
	*/

	static T rootOfSquareSum(double* E1, double* E2, int n){
		T sum = 0;
		for (int i = 0; i < n; i++){
			sum += (E1[i] - E2[i])*(E1[i] - E2[i]);
		}
		return sqrt(sum);
	}

};

#endif