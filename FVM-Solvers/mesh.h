#ifndef MESH_H_H
#define MESH_H_H


#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>
#include "tool.h"
#include "serviceLocator.h"
using namespace std;



class Mesh :public interface_t{
public:

	string meshFile;
	//base
	vector<double> X;
	vector<int> FP;
	vector<int> CP;
	vector<int> IF;//face begin index
	int innerFaceBegin;

	int faceNum;
	int innerFaceNum;
	int cellNum;
	int boundaryFaceNum;
	int nodeNum;
	int boundaryPatchNum;

	vector<int> IC;//cell face begin

	//topology
	vector<int> FC;
	vector<int> CF;

	//geometry
	vector<double> F_a;
	vector<double> F_d;
	vector<double> C_c;
	vector<double> F_c;
	vector<double> F_n;
	vector<double> C_v;
	vector<double> F_e;




	//compress sparse row
	vector<int> IA;
	vector<int> JA;
	vector<int> ON;


	virtual void buildMeshData() = 0;


	//structure
	int Nx, Ny;

};








class GmeshTriangleMesh :public Mesh, public member_t<GmeshTriangleMesh>
{
public:
	//base
	vector<int> F_sort;

	virtual void buildMeshData(){
		ReadFromGshFile(meshFile);
		sortBoundaryFace();
		buildFace();
		generateCompressedSparseRow();
		generateGeo();
	}

	void ReadFromGshFile(string fileName){
		ifstream input(fileName);
		string line;
		while (true){
			getline(input, line);
			if (input.eof())
			{
				break;
			}

			if (line == "$PhysicalNames"){
				boundaryPatchNum = 0;
				int physicalNum;
				input >> physicalNum;
				for (int p = 0; p < physicalNum; p++){
					int d;
					int numberId;
					string name;
					input >> d >> numberId >> name;
					if (d == 1){
						boundaryPatchNum++;
					}
				}
			}

			if (line == "$Nodes"){
				input >> nodeNum;
				X.resize(nodeNum * 2);
				for (int n = 0; n < nodeNum; n++){
					int orderNum;
					double x;
					double y;
					double z;
					input >> orderNum >> x >> y >> z;
					X[n * 2 + 0] = x;
					X[n * 2 + 1] = y;
				}
			}


			if (line == "$Elements"){
				IF.resize(boundaryPatchNum + 2);
				int elementNum;
				input >> elementNum;
				for (int e = 0; e < elementNum; e++){
					int orderNum, elementType, reserveWord, physicalNumId, baseUnitId;
					input >> orderNum >> elementType >> reserveWord >> physicalNumId >> baseUnitId;
					if (elementType == 1){
						int p0, p1;
						input >> p0 >> p1;
						p0--;
						p1--;
						FP.insert(FP.begin() + IF[physicalNumId] * 2 + 0, p0);
						FP.insert(FP.begin() + IF[physicalNumId] * 2 + 1, p1);
						for (int i = physicalNumId; i < boundaryPatchNum + 1; i++){
							IF[i]++;
						}
					}
					else if (elementType == 2){
						int p0, p1, p2;
						input >> p0 >> p1 >> p2;
						p0--;
						p1--;
						p2--;
						CP.push_back(p0);
						CP.push_back(p1);
						CP.push_back(p2);
					}
				}
			}
		}


		boundaryFaceNum = FP.size() / 2;
		cellNum = CP.size() / 3;
	}

	void sortBoundaryFace(){
		int* _RE = new int[boundaryFaceNum * 3];
		for (int f = 0; f < boundaryFaceNum; f++){
			int v0 = FP[f * 2 + 0];
			int v1 = FP[f * 2 + 1];
			int v_min = v0 < v1 ? v0 : v1;
			int v_max = v0>v1 ? v0 : v1;
			_RE[f * 3 + 0] = v_min;
			_RE[f * 3 + 1] = v_max;
			_RE[f * 3 + 2] = f;
		}
		Quick_Sort_RE(_RE, 0, boundaryFaceNum - 1);
		F_sort.resize(boundaryFaceNum);
		for (int f = 0; f < boundaryFaceNum; f++){
			F_sort[f] = _RE[f * 3 + 2];
		}
		delete[]_RE;
	}
	void buildFace(){
		int *_RE = new int[cellNum * 9];	//每条边的相邻面片(v0 v1 c)，其中v0<v1，若有多个相邻面片则重复存储该边。
		for (int c = 0; c<cellNum; c++)
		{
			int v0 = CP[c * 3 + 0];
			int v1 = CP[c * 3 + 1];
			int v2 = CP[c * 3 + 2];

			if (v0<v1)	{ _RE[c * 9 + 0] = v0; _RE[c * 9 + 1] = v1; _RE[c * 9 + 2] = c; }
			else		{ _RE[c * 9 + 0] = v1; _RE[c * 9 + 1] = v0; _RE[c * 9 + 2] = c; }
			if (v1<v2)	{ _RE[c * 9 + 3] = v1; _RE[c * 9 + 4] = v2; _RE[c * 9 + 5] = c; }
			else		{ _RE[c * 9 + 3] = v2; _RE[c * 9 + 4] = v1; _RE[c * 9 + 5] = c; }
			if (v2<v0)	{ _RE[c * 9 + 6] = v2; _RE[c * 9 + 7] = v0; _RE[c * 9 + 8] = c; }
			else		{ _RE[c * 9 + 6] = v0; _RE[c * 9 + 7] = v2; _RE[c * 9 + 8] = c; }
		}

		//Quicksort
		Quick_Sort_RE(_RE, 0, cellNum * 3 - 1);

		int inner_f = 0;
		int boundary_f = 0;
		FC.resize(boundaryFaceNum * 2);
		CF.resize(cellNum * 3);
		for (int i = 0; i<cellNum * 3; i++)
		{
			//printf("edge: %d, %d\n", _RE[i*3+0], _RE[i*3+1]);

			if (i != 3 * cellNum - 1 && _RE[i * 3] == _RE[(i + 1) * 3] && _RE[i * 3 + 1] == _RE[(i + 1) * 3 + 1]) //若当前边和后面一条边相同
			{
				int c0 = _RE[(i + 1) * 3 + 2];
				int c1 = _RE[i * 3 + 2];
				int owner = c0 < c1 ? c0 : c1;//owner 取索引较小的
				int neighbor = c0 > c1 ? c0 : c1;//neighbor 取索引较大的
				//Add the cell to FC
				FC.push_back(owner);
				FC.push_back(neighbor);
				FP.push_back(_RE[i * 3 + 0]);
				FP.push_back(_RE[i * 3 + 1]);
				//Add the face to CF
				CF_on(i, inner_f, _RE);
				CF_on(i + 1, inner_f, _RE);

				inner_f++;
				i++;
			}

			else
			{
				//Add the edge to ET
				FC[F_sort[boundary_f] * 2 + 0] = _RE[i * 3 + 2];
				FC[F_sort[boundary_f] * 2 + 1] = -1;
				//Add the edge to  TE
				int v0 = CP[_RE[i * 3 + 2] * 3 + 0];
				int v1 = CP[_RE[i * 3 + 2] * 3 + 1];
				int v2 = CP[_RE[i * 3 + 2] * 3 + 2];
				if (v0 == _RE[i * 3 + 0] && v1 == _RE[i * 3 + 1] || v1 == _RE[i * 3 + 0] && v0 == _RE[i * 3 + 1]){
					int local_f = _RE[i * 3 + 2] * 3 + 0;
					CF[local_f] = F_sort[boundary_f];
				}
				if (v1 == _RE[i * 3 + 0] && v2 == _RE[i * 3 + 1] || v2 == _RE[i * 3 + 0] && v1 == _RE[i * 3 + 1]){
					int local_f = _RE[i * 3 + 2] * 3 + 1;
					CF[local_f] = F_sort[boundary_f];
				}
				if (v2 == _RE[i * 3 + 0] && v0 == _RE[i * 3 + 1] || v0 == _RE[i * 3 + 0] && v2 == _RE[i * 3 + 1]){
					int local_f = _RE[i * 3 + 2] * 3 + 2;
					CF[local_f] = F_sort[boundary_f];
				}
				boundary_f++;
			}

		}
		delete[]_RE;
		faceNum = FP.size() / 2;
		innerFaceBegin = boundaryPatchNum;
		innerFaceNum = faceNum - boundaryFaceNum;
		IF[IF.size() - 1] = faceNum;

	}

	void CF_on(int i, int inner_f, int* _RE){
		int v0 = CP[_RE[i * 3 + 2] * 3 + 0];
		int v1 = CP[_RE[i * 3 + 2] * 3 + 1];
		int v2 = CP[_RE[i * 3 + 2] * 3 + 2];
		if (v0 == _RE[i * 3 + 0] && v1 == _RE[i * 3 + 1] || v1 == _RE[i * 3 + 0] && v0 == _RE[i * 3 + 1]){
			int local_f = _RE[i * 3 + 2] * 3 + 0;
			CF[local_f] = boundaryFaceNum + inner_f;
		}
		if (v1 == _RE[i * 3 + 0] && v2 == _RE[i * 3 + 1] || v2 == _RE[i * 3 + 0] && v1 == _RE[i * 3 + 1]){
			int local_f = _RE[i * 3 + 2] * 3 + 1;
			CF[local_f] = boundaryFaceNum + inner_f;
		}
		if (v2 == _RE[i * 3 + 0] && v0 == _RE[i * 3 + 1] || v0 == _RE[i * 3 + 0] && v2 == _RE[i * 3 + 1]){
			int local_f = _RE[i * 3 + 2] * 3 + 2;
			CF[local_f] = boundaryFaceNum + inner_f;
		}
	}


	void Quick_Sort_RE(int a[], int l, int r)
	{
		if (l<r)
		{
			int j = Quick_Sort_Partition_RE(a, l, r);

			Quick_Sort_RE(a, l, j - 1);
			Quick_Sort_RE(a, j + 1, r);
		}
	}

	int Quick_Sort_Partition_RE(int a[], int l, int r)
	{
		int pivot[3], i, j, c[3];
		pivot[0] = a[l * 3 + 0];
		pivot[1] = a[l * 3 + 1];
		pivot[2] = a[l * 3 + 2];
		i = l; j = r + 1;
		while (1)
		{
			do ++i; while ((a[i * 3]<pivot[0] || a[i * 3] == pivot[0] && a[i * 3 + 1] <= pivot[1]) && i <= r);
			do --j; while (a[j * 3]>pivot[0] || a[j * 3] == pivot[0] && a[j * 3 + 1]> pivot[1]);
			if (i >= j) break;
			//Swap i and j			
			c[0] = a[i * 3 + 0];
			c[1] = a[i * 3 + 1];
			c[2] = a[i * 3 + 2];
			a[i * 3 + 0] = a[j * 3 + 0];
			a[i * 3 + 1] = a[j * 3 + 1];
			a[i * 3 + 2] = a[j * 3 + 2];
			a[j * 3 + 0] = c[0];
			a[j * 3 + 1] = c[1];
			a[j * 3 + 2] = c[2];
		}
		//Swap l and j
		c[0] = a[l * 3 + 0];
		c[1] = a[l * 3 + 1];
		c[2] = a[l * 3 + 2];
		a[l * 3 + 0] = a[j * 3 + 0];
		a[l * 3 + 1] = a[j * 3 + 1];
		a[l * 3 + 2] = a[j * 3 + 2];
		a[j * 3 + 0] = c[0];
		a[j * 3 + 1] = c[1];
		a[j * 3 + 2] = c[2];
		return j;
	}




	void generateCompressedSparseRow(){
		IA.resize(cellNum + 1);

		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			IA[owner + 1]++;
			IA[neighbor + 1]++;
		}
		for (int i = 1; i<cellNum + 1; i++){
			IA[i]++;
			IA[i] += IA[i - 1];
		}
		JA.resize(IA[cellNum]);
		vector<int> count(cellNum);
		ON.resize(innerFaceNum * 2);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			JA[IA[owner]] = owner;
			count[owner]++;
			ON[f * 2 + 0] = IA[owner] + count[owner];
			JA[ON[f * 2 + 0]] = neighbor;

			JA[IA[neighbor]] = neighbor;
			count[neighbor]++;
			ON[f * 2 + 1] = IA[neighbor] + count[neighbor];
			JA[ON[f * 2 + 1]] = owner;
		}
	}


	void generateGeo(){
		C_c.resize(cellNum * 2);
		C_v.resize(cellNum);
		//special for triangle shape
		for (int c = 0; c<cellNum; c++){
			int p0 = CP[c * 3 + 0];
			int p1 = CP[c * 3 + 1];
			int p2 = CP[c * 3 + 2];
			double x0 = X[p0 * 2 + 0];
			double x1 = X[p1 * 2 + 0];
			double x2 = X[p2 * 2 + 0];
			double xc = 0.25*(x0 + x1 + x2);
			C_c[c * 2 + 0] = xc;
			double y0 = X[p0 * 2 + 1];
			double y1 = X[p1 * 2 + 1];
			double y2 = X[p2 * 2 + 1];
			double yc = 0.25*(y0 + y1 + y2);
			C_c[c * 2 + 1] = yc;

			C_v[c] = VectorMath<double>::triangleArea(x0 - x2, y0 - y2, x1 - x2, y1 - y2);

		}
		//Face
		F_c.resize(faceNum * 2);
		F_a.resize(faceNum);
		F_n.resize(faceNum * 2);


		for (int f = 0; f<faceNum; f++){
			int p1 = FP[f * 2 + 0];
			int p2 = FP[f * 2 + 1];
			double x0 = X[p1 * 2 + 0];
			double x1 = X[p2 * 2 + 0];
			F_c[f * 2 + 0] = 0.5*(x0 + x1);
			double y0 = X[p1 * 2 + 1];
			double y1 = X[p2 * 2 + 1];
			F_c[f * 2 + 1] = 0.5*(y0 + y1);
			F_a[f] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));

			double* orth = VectorMath<double>::coorthogonal(x0 - x1, y0 - y1);
			int ownner = FC[f * 2 + 0];
			double cx0 = C_c[ownner * 2 + 0];
			double cy0 = C_c[ownner * 2 + 1];
			double cx1 = F_c[f * 2 + 0];
			double cy1 = F_c[f * 2 + 1];
			double* norm = VectorMath<double>::normalize(orth[0], orth[1]);
			if (VectorMath<double>::dot(cx1 - cx0, cy1 - cy0, orth[0], orth[1])<0){
				norm = VectorMath<double>::reverse(norm[0], norm[1]);
			}
			F_n[f * 2 + 0] = norm[0];
			F_n[f * 2 + 1] = norm[1];
			delete[] norm;
			delete[] orth;
		}

		//F_d
		F_d.resize(faceNum);
		F_e.resize(faceNum * 2);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int ownner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			double x0 = C_c[ownner * 2 + 0];
			double y0 = C_c[ownner * 2 + 1];
			double x1 = C_c[neighbor * 2 + 0];
			double y1 = C_c[neighbor * 2 + 1];
			F_d[innerF] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			F_e[innerF * 2 + 0] = (x1 - x0) / F_d[innerF];
			F_e[innerF * 2 + 1] = (y1 - y0) / F_d[innerF];
		}
		for (int p = 0; p < boundaryPatchNum; p++){
			for (int bf = IF[p]; bf < IF[p + 1]; bf++){
				int ownner = FC[bf * 2 + 0];
				double x0 = C_c[ownner * 2 + 0];
				double y0 = C_c[ownner * 2 + 1];
				double x1 = F_c[bf * 2 + 0];
				double y1 = F_c[bf * 2 + 1];
				F_d[bf] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
				F_e[bf * 2 + 0] = (x1 - x0) / F_d[bf];
				F_e[bf * 2 + 1] = (y1 - y0) / F_d[bf];
			}
		}
	}




};






//to be done!

class VersatileGmeshMesh :public Mesh, public member_t<VersatileGmeshMesh>
{
public:
	//base
	vector<int> F_sort;



	virtual void buildMeshData(){
		ReadFromGshFile(meshFile);
		sortBoundaryFace();
		buildFace();
		generateCompressedSparseRow();
		generateGeo();
	}

	void ReadFromGshFile(string fileName){
		ifstream input(fileName);
		string line;
		while (true){
			getline(input, line);
			if (input.eof())
			{
				break;
			}

			if (line == "$PhysicalNames"){
				boundaryPatchNum = 0;
				int physicalNum;
				input >> physicalNum;
				for (int p = 0; p < physicalNum; p++){
					int d;
					int numberId;
					string name;
					input >> d >> numberId >> name;
					if (d == 1){
						boundaryPatchNum++;
					}
				}
			}

			if (line == "$Nodes"){
				input >> nodeNum;
				X.resize(nodeNum * 2);
				for (int n = 0; n < nodeNum; n++){
					int orderNum;
					double x;
					double y;
					double z;
					input >> orderNum >> x >> y >> z;
					X[n * 2 + 0] = x;
					X[n * 2 + 1] = y;
				}
			}


			if (line == "$Elements"){
				IC.push_back(0);
				IF.resize(boundaryPatchNum + 2);
				int elementNum;
				input >> elementNum;
				cellNum = 0;
				for (int e = 0; e < elementNum; e++){
					int orderNum, elementType, reserveWord, physicalNumId, baseUnitId;
					input >> orderNum >> elementType >> reserveWord >> physicalNumId >> baseUnitId;
					if (elementType == 1){

						int p0, p1;
						input >> p0 >> p1;
						p0--;
						p1--;
						FP.insert(FP.begin() + IF[physicalNumId] * 2 + 0, p0);
						FP.insert(FP.begin() + IF[physicalNumId] * 2 + 1, p1);
						for (int i = physicalNumId; i < boundaryPatchNum + 1; i++){
							IF[i]++;
						}
					}
					else if (elementType == 2){
						IC.push_back(3);
						int p0, p1, p2;
						input >> p0 >> p1 >> p2;
						p0--;
						p1--;
						p2--;
						CP.push_back(p0);
						CP.push_back(p1);
						CP.push_back(p2);
						cellNum++;
					}
					else if (elementType = 3){
						IC.push_back(4);
						int p0, p1, p2, p3;
						input >> p0 >> p1 >> p2 >> p3;
						p0--;
						p1--;
						p2--;
						p3--;
						CP.push_back(p0);
						CP.push_back(p1);
						CP.push_back(p2);
						CP.push_back(p3);
						cellNum++;
					}
				}
			}
		}
		for (int i = 1; i < IC.size(); i++){
			IC[i] += IC[i - 1];
		}
		boundaryFaceNum = FP.size() / 2;
		innerFaceBegin = boundaryPatchNum;
	}

	void sortBoundaryFace(){
		int* _RE = new int[boundaryFaceNum * 3];
		for (int f = 0; f < boundaryFaceNum; f++){
			int v0 = FP[f * 2 + 0];
			int v1 = FP[f * 2 + 1];
			int v_min = v0 < v1 ? v0 : v1;
			int v_max = v0>v1 ? v0 : v1;
			_RE[f * 3 + 0] = v_min;
			_RE[f * 3 + 1] = v_max;
			_RE[f * 3 + 2] = f;
		}
		Quick_Sort_RE(_RE, 0, boundaryFaceNum - 1);
		F_sort.resize(boundaryFaceNum);
		for (int f = 0; f < boundaryFaceNum; f++){
			F_sort[f] = _RE[f * 3 + 2];
		}
		delete[]_RE;
	}

	void buildFace(){
		int *_RE;//每条边的相邻面片(v0 v1 c)，其中v0<v1，若有多个相邻面片则重复存储该边。
		vector<int> _RE_v;
		for (int c = 0; c<cellNum; c++)
		{
			int cellfaceNum = IC[c + 1] - IC[c];
			int* v = new int[cellfaceNum];
			for (int i = 0; i < cellfaceNum; i++){
				for (int j = i + 1; j < cellfaceNum; j++){
					v[i] = CP[IC[c] + i];
					v[j] = CP[IC[c] + j];
					if (v[i] < v[j]){
						_RE_v.push_back(v[i]);
						_RE_v.push_back(v[j]);
					}
					else{
						_RE_v.push_back(v[j]);
						_RE_v.push_back(v[i]);
					}
					_RE_v.push_back(c);
				}
			}
			delete[] v;
		}
		_RE = _RE_v.data();
		//Quicksort
		Quick_Sort_RE(_RE, 0, _RE_v.size() / 3 - 1);

		int inner_f = 0;
		int boundary_f = 0;
		FC.resize(boundaryFaceNum * 2);
		CF.resize(IC[IC.size() - 1]);
		vector<int> CF_edge_flag(cellNum, 0);
		for (int i = 0; i<_RE_v.size() / 3 - 1; i++)
		{
			//printf("edge: %d, %d\n", _RE[i*3+0], _RE[i*3+1]);

			if (i != _RE_v.size() / 3 - 1 && _RE[i * 3] == _RE[(i + 1) * 3] && _RE[i * 3 + 1] == _RE[(i + 1) * 3 + 1]) //若当前边和后面一条边相同
			{
				int c0 = _RE[(i + 1) * 3 + 2];
				int c1 = _RE[i * 3 + 2];
				int owner = c0 < c1 ? c0 : c1;//owner 取索引较小的
				int neighbor = c0 > c1 ? c0 : c1;//neighbor 取索引较大的
				//Add the cell to FC

				FC.push_back(owner);
				FC.push_back(neighbor);
				FP.push_back(_RE[i * 3 + 0]);
				FP.push_back(_RE[i * 3 + 1]);
				//Add the face to CF
				CF_on(_RE, i, CF_edge_flag.data(), boundaryFaceNum + inner_f);
				CF_on(_RE, i + 1, CF_edge_flag.data(), boundaryFaceNum + inner_f);
				//CF_on(i, inner_f, _RE);
				//CF_on(i + 1, inner_f, _RE);

				inner_f++;
				i++;
			}

			else if (boundary_f<F_sort.size() && FP[F_sort[boundary_f] * 2 + 0] == _RE[i * 3 + 0] && FP[F_sort[boundary_f] * 2 + 1] == _RE[i * 3 + 1]
				|| boundary_f<F_sort.size() && FP[F_sort[boundary_f] * 2 + 1] == _RE[i * 3 + 0] && FP[F_sort[boundary_f] * 2 + 0] == _RE[i * 3 + 1]
				)
			{
				//Add the edge to FC
				FC[F_sort[boundary_f] * 2 + 0] = _RE[i * 3 + 2];
				FC[F_sort[boundary_f] * 2 + 1] = -1;
				//Add the edge to  CF
				CF_on(_RE, i, CF_edge_flag.data(), F_sort[boundary_f]);
				boundary_f++;
			}

		}
		faceNum = FP.size() / 2;
		innerFaceNum = faceNum - boundaryFaceNum;
		IF[IF.size() - 1] = faceNum;

	}

	void CF_on(int* _RE, int e, int* CF_edge_flag, int face){
		int c = _RE[e * 3 + 2];
		int cellfaceNum = IC[c + 1] - IC[c];
		for (int i = 0; i < cellfaceNum; i++){
			for (int j = i + 1; j < cellfaceNum; j++){
				int cp0 = CP[IC[c] + i];
				int cp1 = CP[IC[c] + j];
				int v0 = _RE[e * 3 + 0];
				int v1 = _RE[e * 3 + 1];
				if (cp0 == v0 && cp1 == v1 || cp1 == v0 && cp0 == v1){
					CF[IC[c] + CF_edge_flag[c]] = face;//boundaryFaceNum + inner_f;
					CF_edge_flag[c]++;
					goto breakLoop;
				}
			}
		}
	breakLoop:return;
	}

	void Quick_Sort_RE(int a[], int l, int r)
	{
		if (l<r)
		{
			int j = Quick_Sort_Partition_RE(a, l, r);

			Quick_Sort_RE(a, l, j - 1);
			Quick_Sort_RE(a, j + 1, r);
		}
	}

	int Quick_Sort_Partition_RE(int a[], int l, int r)
	{
		int pivot[3], i, j, c[3];
		pivot[0] = a[l * 3 + 0];
		pivot[1] = a[l * 3 + 1];
		pivot[2] = a[l * 3 + 2];
		i = l; j = r + 1;
		while (1)
		{
			do ++i; while ((a[i * 3]<pivot[0] || a[i * 3] == pivot[0] && a[i * 3 + 1] <= pivot[1]) && i <= r);
			do --j; while (a[j * 3]>pivot[0] || a[j * 3] == pivot[0] && a[j * 3 + 1]> pivot[1]);
			if (i >= j) break;
			//Swap i and j			
			c[0] = a[i * 3 + 0];
			c[1] = a[i * 3 + 1];
			c[2] = a[i * 3 + 2];
			a[i * 3 + 0] = a[j * 3 + 0];
			a[i * 3 + 1] = a[j * 3 + 1];
			a[i * 3 + 2] = a[j * 3 + 2];
			a[j * 3 + 0] = c[0];
			a[j * 3 + 1] = c[1];
			a[j * 3 + 2] = c[2];
		}
		//Swap l and j
		c[0] = a[l * 3 + 0];
		c[1] = a[l * 3 + 1];
		c[2] = a[l * 3 + 2];
		a[l * 3 + 0] = a[j * 3 + 0];
		a[l * 3 + 1] = a[j * 3 + 1];
		a[l * 3 + 2] = a[j * 3 + 2];
		a[j * 3 + 0] = c[0];
		a[j * 3 + 1] = c[1];
		a[j * 3 + 2] = c[2];
		return j;
	}




	void generateCompressedSparseRow(){
		IA.resize(cellNum + 1);

		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			IA[owner + 1]++;
			IA[neighbor + 1]++;
		}
		for (int i = 1; i<cellNum + 1; i++){
			IA[i]++;
			IA[i] += IA[i - 1];
		}
		JA.resize(IA[cellNum]);
		vector<int> count(cellNum);
		ON.resize(innerFaceNum * 2);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			JA[IA[owner]] = owner;
			count[owner]++;
			ON[f * 2 + 0] = IA[owner] + count[owner];
			JA[ON[f * 2 + 0]] = neighbor;

			JA[IA[neighbor]] = neighbor;
			count[neighbor]++;
			ON[f * 2 + 1] = IA[neighbor] + count[neighbor];
			JA[ON[f * 2 + 1]] = owner;
		}
	}


	void generateGeo(){
		C_c.resize(cellNum * 2);
		C_v.resize(cellNum);
		//volumn
		for (int c = 0; c < cellNum; c++){
			double xc = 0;
			double yc = 0;
			for (int p = IC[c]; p < IC[c + 1]; p++){
				xc += X[CP[p] * 2 + 0];
				yc += X[CP[p] * 2 + 1];
			}
			xc /= (IC[c + 1] - IC[c]);
			yc /= (IC[c + 1] - IC[c]);
			C_c[c * 2 + 0] = xc;
			C_c[c * 2 + 1] = yc;
			double vol = 0;
			for (int f = IC[c]; f < IC[c + 1]; f++){

				int p0 = FP[CF[f] * 2 + 0];
				int p1 = FP[CF[f] * 2 + 1];
				double x0 = X[p0 * 2 + 0];
				double x1 = X[p1 * 2 + 0];
				double y0 = X[p0 * 2 + 1];
				double y1 = X[p1 * 2 + 1];
				vol += VectorMath<double>::triangleArea(x0 - xc, y0 - yc, x1 - xc, y1 - yc);
			}
			C_v[c] = vol;
		}
		//Face
		F_c.resize(faceNum * 2);
		F_a.resize(faceNum);
		F_n.resize(faceNum * 2);


		for (int f = 0; f<faceNum; f++){
			int p1 = FP[f * 2 + 0];
			int p2 = FP[f * 2 + 1];
			double x0 = X[p1 * 2 + 0];
			double x1 = X[p2 * 2 + 0];
			F_c[f * 2 + 0] = 0.5*(x0 + x1);
			double y0 = X[p1 * 2 + 1];
			double y1 = X[p2 * 2 + 1];
			F_c[f * 2 + 1] = 0.5*(y0 + y1);
			F_a[f] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));

			double* orth = VectorMath<double>::coorthogonal(x0 - x1, y0 - y1);
			int ownner = FC[f * 2 + 0];
			double cx0 = C_c[ownner * 2 + 0];
			double cy0 = C_c[ownner * 2 + 1];
			double cx1 = F_c[f * 2 + 0];
			double cy1 = F_c[f * 2 + 1];
			double* norm = VectorMath<double>::normalize(orth[0], orth[1]);
			if (VectorMath<double>::dot(cx1 - cx0, cy1 - cy0, orth[0], orth[1])<0){
				norm = VectorMath<double>::reverse(norm[0], norm[1]);
			}
			F_n[f * 2 + 0] = norm[0];
			F_n[f * 2 + 1] = norm[1];
			delete[] norm;
			delete[] orth;
		}

		//F_d
		F_d.resize(faceNum);
		F_e.resize(faceNum * 2);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = IF[innerFaceBegin] + f;
			int ownner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			double x0 = C_c[ownner * 2 + 0];
			double y0 = C_c[ownner * 2 + 1];
			double x1 = C_c[neighbor * 2 + 0];
			double y1 = C_c[neighbor * 2 + 1];
			F_d[innerF] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			F_e[innerF * 2 + 0] = (x1 - x0) / F_d[innerF];
			F_e[innerF * 2 + 1] = (y1 - y0) / F_d[innerF];
		}
		for (int p = 0; p < boundaryPatchNum; p++){
			for (int bf = IF[p]; bf < IF[p + 1]; bf++){
				int ownner = FC[bf * 2 + 0];
				double x0 = C_c[ownner * 2 + 0];
				double y0 = C_c[ownner * 2 + 1];
				double x1 = F_c[bf * 2 + 0];
				double y1 = F_c[bf * 2 + 1];
				F_d[bf] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
				F_e[bf * 2 + 0] = (x1 - x0) / F_d[bf];
				F_e[bf * 2 + 1] = (y1 - y0) / F_d[bf];
			}
		}
	}
};



#endif