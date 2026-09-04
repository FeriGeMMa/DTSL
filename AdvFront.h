#pragma once
#pragma warning(disable : 4996)
#include"vector"


typedef struct DoubleLinkedList1
{
	double key, y;
	long TriangleIndex; // pointer to triangle
	long EntryIndex;
	long VertexNo;      // vertexIndex

	DoubleLinkedList1* RPointer;
	DoubleLinkedList1* LPointer;

	DoubleLinkedList1* Next;   // for making the stack
} DoubleLinkedList;




class AdvFront
{
private:
	double XMax, XMin, Size;
	double RoundingTolerance, SmallTriangleTolerance;
	int NoOfCells;
	std::vector<DoubleLinkedList*> SkipTable;

	DoubleLinkedList* Left;
	DoubleLinkedList* Right;
	DoubleLinkedList* CurrentLeft;
	DoubleLinkedList* CurrentRight;

	DoubleLinkedList1* Pool;

	int ReturnEntryIndex(double x);

	DoubleLinkedList* ReturnFromPool();
	void InsertIntoPool(DoubleLinkedList* w);


public:

	void Init(int Fraction, double xmin, double ymin, double xmax, double rt);

	DoubleLinkedList* ReturnLeft() { return Left; }
	DoubleLinkedList* ReturnRight() { return Right; }

	DoubleLinkedList* ReturnCurrentLeft() { return CurrentLeft; }
	DoubleLinkedList* ReturnCurrentRight() { return CurrentRight; }


	DoubleLinkedList* Insert(double x, double y, long vertex, long TIndex);
	void InsertRightFrom(DoubleLinkedList* RefPointer, double x, double y, int v, int ti);

	int Projection(double x, double y);
	int UpdateTIndex(double x, long i, long j, long ti);
	void DeleteElement(DoubleLinkedList* Pointer, int TIndex);

	void ClearAF();
	int ClearAFInBetween();


	// aux
	void PrintAF(const char* fname);
	void PrintAFCoordinates(const char* fname);
};

