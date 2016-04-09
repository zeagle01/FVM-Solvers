#ifndef PRINTER_H_H
#define PRINTER_H_H


#include "mesh.h"
#include <fstream>
#include <vector>
#include "dataStructure.h"
#include "serviceLocator.h"
#include <sstream>
using namespace std;
class  A_Printer :public interface_t{

public:
	string cases = "cases ";
	string zoneType;
	virtual void print(string fileName, Mesh* mesh, vector<CellField> phi) = 0;
};

class  Printer :public A_Printer, public member_t<Printer>{
public:
	void printFileHead(ofstream &fout, Mesh* mesh, vector<CellField> phi){

		fout << "variables = \"x\", \"y\"," << " ";
		for (int vn = 0; vn < phi.size(); vn++){
			fout << " \"v" << vn << "\",";
		}
		fout << endl;
		fout << "ZONE I = " << mesh->Ny << " J = " << mesh->Nx << " F = POINT";
		fout << endl;

	}
	void printBoundary(ofstream &fout, Mesh* mesh, vector<CellField> phi){

		for (int f = 0; f < mesh->boundaryFaceNum; f++){
			int boundaryF = mesh->IF[0] + f;
			fout << mesh->F_c[boundaryF * 2 + 0] << " " << mesh->F_c[boundaryF * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn].boundary[f] << " ";
			}
			fout << endl;
		}

	}
	void printInner(ofstream &fout, Mesh* mesh, vector<CellField> phi){
		for (int c = 0; c < mesh->cellNum; ++c){
			fout << mesh->C_c[c * 2 + 0] << " " << mesh->C_c[c * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn].inner[c] << " ";
			}
			fout << endl;
		}
	}
	void print(string fileName, Mesh* mesh, vector<CellField> phi)
	{
		ofstream fout(cases + fileName);

		printFileHead(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};


class  Printer1 :public Printer, public member_t<Printer>{
public:
	void print(string fileName, Mesh* mesh, vector<CellField> phi)
	{
		ofstream fout(cases + fileName);
		printFileHead(fout, mesh, phi);
		printBoundary(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};



class  UnstructurePrinter :public A_Printer, public member_t<UnstructurePrinter>{
public:
	void printFileHead(ofstream &fout, Mesh* mesh, vector<CellField> phi){

		fout << "variables = \"x\", \"y\"," << " ";
		for (int vn = 0; vn < phi.size(); vn++){
			fout << " \"v" << vn << "\",";
		}
		fout << endl;
		fout << "ZONE N = " << mesh->nodeNum << " E = " << mesh->cellNum << " ZONETYPE =fe" << zoneType << endl;;
		fout << "DATAPACKING=BLOCK" << endl;
		stringstream ss;
		ss << "-" << 2 + phi.size();
		string s1 = ss.str();
		string s = phi.size() == 1 ? ("") :s1;
		fout << "VARLOCATION=([1-2]=NODAL,[3" << s << "]=CELLCENTERED)" << endl;
		fout << endl;

	}
	void printBoundary(ofstream &fout, Mesh* mesh, vector<CellField> phi){

		for (int f = 0; f < mesh->boundaryFaceNum; f++){
			int boundaryF = mesh->IF[0] + f;
			fout << mesh->F_c[boundaryF * 2 + 0] << " " << mesh->F_c[boundaryF * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn].boundary[f] << " ";
			}
			fout << endl;
		}

	}
	void printInner(ofstream &fout, Mesh* mesh, vector<CellField> phi){
		for (int n = 0; n < mesh->nodeNum; ++n){
			fout << mesh->X[n * 2 + 0] << endl;

		}
		fout << endl;
		for (int n = 0; n < mesh->nodeNum; ++n){
			fout << mesh->X[n * 2 + 1] << endl;

		}
		fout << endl;
		for (int vn = 0; vn < phi.size(); vn++){
			for (int c = 0; c < mesh->cellNum; ++c){
				fout << phi[vn].inner[c] << endl;

			}
			fout << endl;
		}
		fout << endl;
		for (int i = 0; i < mesh->cellNum; i++){
			for (int j = mesh->IC[i]; j < mesh->IC[i + 1]; j++){
				fout << mesh->CP[j] + 1 << "	";
			}
			fout << endl;
		}
		fout << endl;
	}
	void print(string fileName, Mesh* mesh, vector<CellField> phi)
	{
		ofstream fout(cases + fileName);

		printFileHead(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};

#endif