#pragma warning(disable : 4996)

#include <stdio.h>
#include "AuxF.h"
#include "AdvFront.h"


int AFSearchCounter = 0;
int AFJumpCounter = 0;

void AdvFront::Init(int fraction, double xmin, double ymin, double xmax, double rt)
{
	RoundingTolerance = rt;
	SmallTriangleTolerance = 10 * rt;

	double dx = 0.0009999999776482582 * (xmax - xmin);

	xmin -= dx;
	xmax += dx;
	ymin -= dx;

	Size = (xmax - xmin) / fraction;
	NoOfCells = int((xmax - xmin) / Size) + 1;

	XMax = xmax;
	XMin = xmin;

	SkipTable.clear();

	SkipTable.resize(NoOfCells);
	for (int i = 0; i < SkipTable.size(); i++)
		SkipTable[i] = nullptr;

	Left = new DoubleLinkedList;
	Left->EntryIndex = 0;
	Left->key = xmin;
	Left->y = ymin;
	Left->VertexNo = 0;
	Left->LPointer = nullptr;
	Left->TriangleIndex = NONE;
	Left->Next = nullptr;

	Right = new DoubleLinkedList;
	Right->EntryIndex = NoOfCells - 1;
	Right->key = xmax;
	Right->y = ymin;
	Right->VertexNo = 1;
	Right->LPointer = Left;
	Right->RPointer = nullptr;
	Right->TriangleIndex = NONE;
	Right->Next = nullptr;
	Left->RPointer = Right;
	

	SkipTable[0] = Left;
	SkipTable[NoOfCells-1] = Right;

	CurrentLeft = Left;
	CurrentRight = Right;

	Pool = nullptr;
}


DoubleLinkedList* AdvFront::Insert(double x, double y, long vertex, long TIndex)
{
	DoubleLinkedList* NewEl = ReturnFromPool();
	int STIndex = ReturnEntryIndex(x);

	NewEl->EntryIndex = STIndex;
	NewEl->key = x;
	NewEl->VertexNo = vertex;
	NewEl->y = y;
	NewEl->TriangleIndex = TIndex;
	NewEl->Next = nullptr;

	CurrentLeft->RPointer = NewEl;
	NewEl->RPointer = CurrentRight;
	CurrentRight->LPointer = NewEl;
	NewEl->LPointer = CurrentLeft;

	CurrentLeft->TriangleIndex = TIndex;

	if (SkipTable[STIndex] == nullptr)
		SkipTable[STIndex] = NewEl;
	else // correct the SkipTable pointer, if necessary

		if ((SkipTable[STIndex]->key > x) || Equal(SkipTable[STIndex]->key, x, RoundingTolerance))
			SkipTable[STIndex] = NewEl;

	return NewEl;
}

DoubleLinkedList* AdvFront::ReturnFromPool()
{
	if (Pool == nullptr)
	{
		DoubleLinkedList* w = new DoubleLinkedList;
		w->Next = nullptr;
		return w;
	}

	DoubleLinkedList* w = Pool;
	Pool = Pool->Next;
	return w;
}

void AdvFront::InsertIntoPool(DoubleLinkedList* w)
{
	if (w != nullptr)
	{
		w->Next = Pool;
		Pool = w;
	}
}

void AdvFront::ClearAF()
{
	DoubleLinkedList* w;
	while (Left != nullptr)
	{
		w = Left;
		Left = Left->RPointer;
		delete w;
	}

	while (Pool != nullptr)
	{
		w = Pool;
		Pool = Pool->Next;
		delete w;
	}
	Pool = nullptr;
}


int AdvFront::ReturnEntryIndex(double x)
{
	int index = (int)((x - XMin) / Size);
	if (index < 0) index = 0;
	return index;

}

int  AdvFront::Projection(double x, double y)
{
	int k = ReturnEntryIndex(x);
	DoubleLinkedList* w;


//	while ((w = SkipTable[k]) == NULL)
//		k++;


	while ((w = SkipTable[k]) == nullptr)
	{
		k++;
		AFJumpCounter++;
	}

//	while (w->key < x)
//		w = w->RPointer;

	while (w->key < x)
	{
		w = w->RPointer;
		AFSearchCounter++;
	}

	CurrentRight = w;
	CurrentLeft = w->LPointer;

	if (TheSamePoints2D(CurrentRight->key, CurrentRight->y, x, y, RoundingTolerance))
		return 0;
	if (TheSamePoints2D(CurrentLeft->key, CurrentLeft->y, x, y, RoundingTolerance))
		return 0;

	if (Pitagora(x, y, CurrentRight->key, CurrentRight->y) <= SmallTriangleTolerance)
		return 0;

	if (Pitagora(x, y, CurrentLeft->key, CurrentLeft->y) <= SmallTriangleTolerance)
		return 0;

	return 1;
}


int AdvFront::UpdateTIndex(double x, long i, long j, long ti)
{
	int EntryIndex = (long)((x - XMin) / Size);
	if (EntryIndex < 0)
		EntryIndex = 0;
	if (EntryIndex >= NoOfCells)
		EntryIndex = NoOfCells - 1;

	while (SkipTable[EntryIndex] == nullptr)
		EntryIndex--;

	DoubleLinkedList* w = SkipTable[EntryIndex];


	while (1)
	{
		if (w->RPointer == nullptr)
			return 0;

		if ((w->VertexNo == i) && (w->RPointer->VertexNo == j))
		{
			w->TriangleIndex = ti;
			return 1;
		}
		if ((w->VertexNo == j) && (w->RPointer->VertexNo == i))
		{
			w->TriangleIndex = ti;
			return 1;
		}

		if ((w->VertexNo == i) && (w->LPointer->VertexNo == j))
		{
			w->TriangleIndex = ti;
			return 1;
		}

		if ((w->VertexNo == j) && (w->LPointer->VertexNo == i))
		{
			w->TriangleIndex = ti;
			return 1;
		}

		if (w->key > x)
			return 0;
		if ((w = w->RPointer) == nullptr)
			return 0;
	}
}


void AdvFront::DeleteElement(DoubleLinkedList* Pointer, int TIndex)
{
	int STIndex = Pointer->EntryIndex;
	if (SkipTable[STIndex] == Pointer) // Skip table point on the element, which will be removed
	{
		DoubleLinkedList* nextR = Pointer->RPointer;
		if (nextR->EntryIndex == STIndex)
			SkipTable[STIndex] = nextR;
		else
			SkipTable[STIndex] = nullptr;
	}

	Pointer->LPointer->RPointer = Pointer->RPointer;
	Pointer->RPointer->LPointer = Pointer->LPointer;
	Pointer->LPointer->TriangleIndex = TIndex;

	InsertIntoPool(Pointer);
}




int AdvFront::ClearAFInBetween()
// Method removes records between the th second first and the second last record
//It dows not affect HashTable pointers as it is used only during the finalization; 
{
	if (Left->RPointer == Right) // AF is empty
		return 0; // 

	DoubleLinkedList* leftw = Left->RPointer;
	
	if ((leftw = leftw->RPointer) == Right->LPointer)  // there is only between the first and the last
		return 0;

	DoubleLinkedList* L = Left->RPointer;
	DoubleLinkedList* R = Right->LPointer;
	DoubleLinkedList* p;
	 

	while ((leftw != nullptr) && (leftw != Right->LPointer))
	{
		p = leftw;
		leftw = leftw->RPointer;
		InsertIntoPool(p);
	}

	L->RPointer = R;
	R->LPointer = L;

	return 1;
}


void AdvFront::InsertRightFrom(DoubleLinkedList* RefPointer, double x, double y, int v, int ti)
// Thiw function does not updates the pointer to the hash table as it is called during the finalisation only
{ 
	DoubleLinkedList* NewEl = ReturnFromPool();

	NewEl->key = x;
	NewEl->y = y;
	NewEl->VertexNo = v;
	NewEl->TriangleIndex = ti;
	NewEl->EntryIndex = ReturnEntryIndex(x);

	NewEl->RPointer = RefPointer->RPointer;
	RefPointer->RPointer->LPointer = NewEl;
	RefPointer->RPointer = NewEl;
	NewEl->LPointer = RefPointer;
}




void AdvFront::PrintAF(const char* fname)
{
	FILE* f = fopen(fname, "wt");
	DoubleLinkedList* w = Left;

	int i = 0;
	fprintf(f, "SkipTable \n");
	for (int i = 0; i < SkipTable.size(); i++)
	{
		if (SkipTable[i] != nullptr)
			fprintf(f, "%d %d \n", i, SkipTable[i]->VertexNo);
		else
			fprintf(f, "%d NULL \n", i);
	}

	fprintf(f, "AF \n");
	while (w != nullptr)
	{
		fprintf(f, "v: %d, ti %d \n", w->VertexNo, w->TriangleIndex);
		w = w->RPointer;
	}
	fclose(f);

}


void AdvFront::PrintAFCoordinates(const char* fname)
{
	FILE* f = fopen(fname, "wt");
	DoubleLinkedList* w = Left;
	while (w != nullptr)
	{
		fprintf(f, "%d   %f %f\n", w->VertexNo, w->key, w->y);
		w = w->RPointer;
	}

	fclose(f);
}

