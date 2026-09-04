#pragma once

#pragma warning(disable : 4996)

#include "AdvFront.h"
#include "AuxF.h"


#include <stdio.h>
#include <vector>

#define LEFT_WALK 1
#define RIGHT_WALK 2


typedef struct VertexStruct1// point to be triangulated
{
	long Index;
	double x, y;

} VertexStruct;


typedef struct bt1
{
	long i, j, k;						// triangle coordinate
	long ij,							// pointers to triangle neighbours
		ki,
		jk;
	bool Eij, Eki, Ejk;					// true if has been legalised already;
	bool AFFlag;                        // true if triangle is part of AF
	bool ArtificialFlag;				// true if triangle is connected to artificial vertex; the flag is set during the finalization
} BaseTriangle;

typedef struct stlStl
{
	double x1, y1;
	double x2, y2;
	double x3, y3;

}STLStyleTriangle;


class MySweepLDT
{
private:
	int NoOfVertices;
	double xmin, ymin;
	double xmax, ymax;
	double RoundingTolerance;
	double AllovedGrooveDepth;
	int NextTriangPosition;
	int SwapCounter;
	int TFillCounter;

	AdvFront AF;

	std::vector<VertexStruct> Vertices;
	std::vector<BaseTriangle> TriangArray;


	void  InsertTriangle(int i, int j, int k, int ij, int jk, int ki);
//	void  SetArtificialVertices();

	// main part
	int InsertTopTriangle(VertexStruct&);
	void SetNeighbour(long T1, long v1, long v2, long TIndex);
	void SetAFFlag(int TIndex);

	// legalization
	void CheckDT(int t1, int t2);
	int PrepareTriangles(int t1, int t2);
	void ShiftVerticesRight(BaseTriangle& t);
	void ReplaceNeighbour(int TIndex, int OldT, int NewNeighbour);
	void  UpdateTriangleIndex(int ti);

	// adding triangles left/right
//	DoubleLinkedList* AddSharpTriangles(DoubleLinkedList* CurrP, int LRF);
	int AddSharpTriangles(DoubleLinkedList* CurrP, int LRF);
	int CloseTriangle(DoubleLinkedList* p1, DoubleLinkedList* p2, DoubleLinkedList* p3, int lrf);
	void SetAllAFFlags(int TIndex);
	int GetSecondTriangleIndex(BaseTriangle& t, int ind);
	int SolveNearToVerticalAF(DoubleLinkedList* Top, int lrf);
	int SetWalkers(DoubleLinkedList* CurrP, DoubleLinkedList*& AFWalker, DoubleLinkedList*& AFWalkerNext, int lrf);


	int  Groove(DoubleLinkedList* AFWalker, DoubleLinkedList* AFWalkerNext, int lrf);


	// finalization
	void Finalize();
	void ConstructBootomAF();
	int FillToConvexHull(bool UpperFlag);
	int TwoVerticesArtificial(int ti);
	int ReturnNeighbouringTriangle(int TIndex, int va, int vb);
	int ThirdTriangleVertex(int TIndex, int v1, int v2);

	int TIsNotArtfTriangle(BaseTriangle& bt);

	// aux
	int PrintTriangles(char* fname);
	int PrintNVertices(int n, char* fname);

public:
	int InitSL(double x1, double y1, double x2, double y2, int LutEntries, const std::vector<double>& X, const std::vector<double>& Y);
	void DoDT();
	int ReturnSwapCounter() { return SwapCounter; };
	int ReturnTFillCounter() { return TFillCounter; };
	void ReturnAllTriangles(std::vector<STLStyleTriangle>& Triangles);
	void ReturnDTTriangles(std::vector<STLStyleTriangle>& Triangles);


};

