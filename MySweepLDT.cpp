
#pragma warning(disable : 4996)
#include <stdio.h>
#include <math.h>

#include "MySweepLDT.h"
#include "AuxF.h"

#include <vector>
#include <algorithm>

#define M_PI 3.14159265358979323846

int MySweepLDT::InitSL(double x1, double y1, double x2, double y2, int LutEntries, const std::vector<double>& X, const std::vector<double>& Y)
{

    NoOfVertices = (int)X.size();
    if (NoOfVertices < 2) return 0;

    double dx = x2 - x1;
    double dy = y2 - y1;

    AllovedGrooveDepth = dy * 0.025; 

    xmin = x1 - dx; ymin = y1 - dy;
    xmax = x2 + dx; ymax = y2 + 0.1 * dy;

    if (dx >= dy)
        RoundingTolerance = 0.000001 * dx;
    else
        RoundingTolerance = 0.000001 * dy;


    if (Equal(dx, 0, RoundingTolerance))
        return 0;

    if (Equal(dy, 0, RoundingTolerance))
        return 0;

    Vertices.clear();
    Vertices.resize(NoOfVertices + 2);

    //    Set artificial vertices
    Vertices[0].x = xmin;
    Vertices[0].y = ymin;
    Vertices[0].Index = 0;

    Vertices[1].x = xmax;
    Vertices[1].y = ymin;
    Vertices[1].Index = 1;

    // add input vertices
    for (long i = 0; i < NoOfVertices; i++)
    {
        Vertices[i + 2].x = X[i];
        Vertices[i + 2].y = Y[i];
        Vertices[i + 2].Index = i;

    }

    NoOfVertices = (int)Vertices.size();


    // sert vertices lexicographically where y is the main sorting direction
    std::sort(Vertices.begin(), Vertices.end(), [](const VertexStruct& a, const VertexStruct& b) {
        if (a.y == b.y)
            return a.x < b.x;  // If y is equal, compare x
        return a.y < b.y;      // Otherwise, compare y
        });



    // assign consequtive indices to the vertices
    for (int i = 2; i < Vertices.size(); i++)
        Vertices[i].Index = i;

    TriangArray.clear();

    TriangArray.resize(3 * NoOfVertices);
    NextTriangPosition = 0;
    InsertTriangle(2, 1, 0, NONE, NONE, NONE);

    AF.Init(LutEntries, x1 - 0.05 * dx, y1 - 0.05 * dy, x2 + 0.05 * dx, RoundingTolerance);
    AF.Insert(Vertices[2].x, Vertices[2].y, Vertices[2].Index, NextTriangPosition);

    NextTriangPosition++;
    SwapCounter = 0;
    TFillCounter = 0;

    return 1;
}



void  MySweepLDT::InsertTriangle(int i, int j, int k, int ij, int jk, int ki)
{
    BaseTriangle t;
    t.i = i;
    t.j = j;
    t.k = k;
    t.ij = ij;
    t.jk = jk;
    t.ki = ki;
    t.Eij = t.Ejk = t.Eki = false;
    t.AFFlag = true;
    t.ArtificialFlag = false;
    TriangArray[NextTriangPosition] = t;
}


// main part ###

void MySweepLDT::DoDT()
{
    for (int i = 3; i < Vertices.size(); i++)

   // for (int i = 3; i < 110; i++)
       //    for (int i = 3; i < 80; i++)
    {
        
 //       if ((i % 100) == 0)
/*        if (i >= 4000)
        {
            FILE* f = fopen("d:\\1000.txt", "wt");
            fprintf(f, "%d\n", i);
            fclose(f);
        }
        
        if (i == 18)
        {
            int stop1 = 0;
        }
        */

  
        if (InsertTopTriangle(Vertices[i]))
        {

            DoubleLinkedList* Left = AF.ReturnCurrentLeft();
            DoubleLinkedList* Right = AF.ReturnCurrentRight();

            int Pair = Left->TriangleIndex;

            DoubleLinkedList* ActualAFPointer = AF.Insert(Vertices[i].x, Vertices[i].y, i, NextTriangPosition);

            CheckDT(Pair, NextTriangPosition);
            NextTriangPosition++;

            while (AddSharpTriangles(ActualAFPointer, RIGHT_WALK));
            while (AddSharpTriangles(ActualAFPointer, LEFT_WALK));

            while (SolveNearToVerticalAF(ActualAFPointer, RIGHT_WALK));
            while (SolveNearToVerticalAF(ActualAFPointer, LEFT_WALK));

            while (AddSharpTriangles(ActualAFPointer, RIGHT_WALK));
            while (AddSharpTriangles(ActualAFPointer, LEFT_WALK));
        }


    }
    //    AF.PrintAFCoordinates("d:\\AFCOORDStart.TXT");
        Finalize();
}



int MySweepLDT::InsertTopTriangle(VertexStruct& v)
{
    // triangles should be orientetd in mathematical positive sense
    if (AF.Projection(v.x, v.y))
    {
        DoubleLinkedList* lt = AF.ReturnCurrentLeft();
        DoubleLinkedList* rt = AF.ReturnCurrentRight();
        int  OldT = lt->TriangleIndex;

        InsertTriangle(v.Index, rt->VertexNo, lt->VertexNo, NONE, OldT, NONE);
        SetNeighbour(OldT, rt->VertexNo, lt->VertexNo, NextTriangPosition);
        SetAFFlag(OldT);
        return 1;
    }
    return 0;
}



void MySweepLDT::SetNeighbour(long T1, long v1, long v2, long TIndex)
// Function sets the new neighbour TIndex of triangle T1 over the
// vertices v1 v2.
{
    if (T1 != NONE)
        if (((TriangArray[T1].i == v1) && (TriangArray[T1].j == v2)) ||
            ((TriangArray[T1].i == v2) && (TriangArray[T1].j == v1)))
            TriangArray[T1].ij = TIndex;
        else
        {
            if (((TriangArray[T1].j == v1) && (TriangArray[T1].k == v2)) ||
                ((TriangArray[T1].j == v2) && (TriangArray[T1].k == v1)))
                TriangArray[T1].jk = TIndex;
            else
                TriangArray[T1].ki = TIndex;
        }
}


void MySweepLDT::SetAFFlag(int TIndex)
{

    if ((TriangArray[TIndex].ij != NONE) && (TriangArray[TIndex].jk != NONE) && (TriangArray[TIndex].ki != NONE))
        TriangArray[TIndex].AFFlag = false;
    else
        TriangArray[TIndex].AFFlag = true;
}



void MySweepLDT::CheckDT(int t1, int t2)
// Triangles t1 and t2 are checked for Delaunay criterium and 
// legalised in neccessary
{

    if ((t1 != NONE) && (t2 != NONE))
        if (t1 != t2)
        {
            long i1, j1, k1, k2;

            if (PrepareTriangles(t1, t2))
            {
                if (PrepareTriangles(t2, t1))
                {
                    i1 = TriangArray[t1].i;
                    j1 = TriangArray[t1].j;
                    k1 = TriangArray[t1].k;
                    k2 = TriangArray[t2].k;

                    double C[3][3];

                    C[0][0] = Vertices[i1].x - Vertices[k2].x;
                    C[0][1] = Vertices[i1].y - Vertices[k2].y;
                    C[0][2] = C[0][0] * C[0][0] + C[0][1] * C[0][1];

                    C[1][0] = Vertices[j1].x - Vertices[k2].x;
                    C[1][1] = Vertices[j1].y - Vertices[k2].y;
                    C[1][2] = C[1][0] * C[1][0] + C[1][1] * C[1][1];

                    C[2][0] = Vertices[k1].x - Vertices[k2].x;
                    C[2][1] = Vertices[k1].y - Vertices[k2].y;
                    C[2][2] = C[2][0] * C[2][0] + C[2][1] * C[2][1];

                    double D = C[0][0] * (C[1][1] * C[2][2] - C[1][2] * C[2][1])
                        + C[0][1] * (C[1][2] * C[2][0] - C[1][0] * C[2][2])
                        + C[0][2] * (C[1][0] * C[2][1] - C[1][1] * C[2][0]);

 //                   if (Equal(D, 0, RoundingTolerance)) goto w;

                    if (D < 0)
                    {
                        SwapCounter++;

                        long n_ij1 = TriangArray[t1].ij;
                        long n_jk1 = TriangArray[t1].jk;
                        long n_ki1 = TriangArray[t1].ki;

                        long n_ij2 = TriangArray[t2].ij;
                        long n_jk2 = TriangArray[t2].jk;
                        long n_ki2 = TriangArray[t2].ki;

                        // new triangles are at the location of the old triangles
                        TriangArray[t1].i = k1;
                        TriangArray[t1].j = i1;
                        TriangArray[t1].k = k2;

                        TriangArray[t2].i = k1;
                        TriangArray[t2].j = k2;
                        TriangArray[t2].k = j1;

                        // update the neighbours
                        TriangArray[t1].ij = n_ki1;
                        TriangArray[t1].jk = n_jk2;
                        TriangArray[t1].ki = t2;

                        TriangArray[t2].ij = t1;
                        TriangArray[t2].jk = n_ki2;
                        TriangArray[t2].ki = n_jk1;

                        if (n_jk2 != NONE)
                            ReplaceNeighbour(n_jk2, t2, t1);

                        if (n_jk1 != NONE)
                            ReplaceNeighbour(n_jk1, t1, t2);

                        SetAFFlag(t1);
                        SetAFFlag(t2);
                        if (TriangArray[t1].AFFlag)
                            UpdateTriangleIndex(t1);

                        if (TriangArray[t2].AFFlag)
                            UpdateTriangleIndex(t2);

                        CheckDT(t1, n_jk2);
                        CheckDT(t2, n_ki2);
                        CheckDT(t1, n_ki1);
                        CheckDT(t2, n_jk1);
                    }
                }
            }
        }
//w:;
}



int MySweepLDT::PrepareTriangles(int t1, int t2)
// this function rotates the triangles t1 and t2 in that way 
// that Delaunay criterion can be applied.
{
    while ((TriangArray[t1].k == TriangArray[t2].i) ||
        (TriangArray[t1].k == TriangArray[t2].j) ||
        (TriangArray[t1].k == TriangArray[t2].k))
        ShiftVerticesRight(TriangArray[t1]);

    if ((TriangArray[t1].i == TriangArray[t2].j) &&
        (TriangArray[t1].j == TriangArray[t2].i))
        return 1;

    ShiftVerticesRight(TriangArray[t2]);

    if ((TriangArray[t1].i == TriangArray[t2].j) &&
        (TriangArray[t1].j == TriangArray[t2].i))
        return 1;

    ShiftVerticesRight(TriangArray[t2]);
    if ((TriangArray[t1].i == TriangArray[t2].j) &&
        (TriangArray[t1].j == TriangArray[t2].i))
        return 1;
    return 0;
}

void MySweepLDT::ShiftVerticesRight(BaseTriangle& t)
{
    long a;

    a = t.k;
    t.k = t.j;
    t.j = t.i;
    t.i = a;

    a = t.ki;
    t.ki = t.jk;
    t.jk = t.ij;
    t.ij = a;
}



void MySweepLDT::ReplaceNeighbour(int TIndex, int OldT, int NewNeighbour)
{
    if (TriangArray[TIndex].ij == OldT)
        TriangArray[TIndex].ij = NewNeighbour;
    else
    {
        if (TriangArray[TIndex].jk == OldT)
            TriangArray[TIndex].jk = NewNeighbour;
        else
        {
            if (TriangArray[TIndex].ki == OldT)
                TriangArray[TIndex].ki = NewNeighbour;
        }
    }
}



void  MySweepLDT::UpdateTriangleIndex(int ti)
{
    int DoneFlag = 0;

    if (TriangArray[ti].ij == NONE)
    {
        DoneFlag = AF.UpdateTIndex(Vertices[TriangArray[ti].i].x,
            TriangArray[ti].i,
            TriangArray[ti].j, ti);
    }

    if (TriangArray[ti].jk == NONE)
    {
        DoneFlag = AF.UpdateTIndex(Vertices[TriangArray[ti].j].x,
            TriangArray[ti].j,
            TriangArray[ti].k, ti);
    }

    if (TriangArray[ti].ki == NONE)
    {
        DoneFlag = AF.UpdateTIndex(Vertices[TriangArray[ti].k].x,
            TriangArray[ti].k,
            TriangArray[ti].i, ti);
    }
}

int MySweepLDT::SetWalkers(DoubleLinkedList* CurrP, DoubleLinkedList*& AFWalker, DoubleLinkedList*& AFWalkerNext, int lrf)
{
    if (lrf == LEFT_WALK)
    {
        AFWalker = CurrP->LPointer;
        if (AFWalker == nullptr)
            return 0;

        if (AFWalker == AF.ReturnLeft())
            return 0;
        AFWalkerNext = AFWalker->LPointer;
        if (AFWalkerNext == AF.ReturnLeft())
            return 0;
    }
    else
    {
        AFWalker = CurrP->RPointer;
        if (AFWalker == nullptr)
            return 0;

        if (AFWalker == AF.ReturnRight())
            return 0;
        AFWalkerNext = AFWalker->RPointer;
        if (AFWalkerNext == AF.ReturnRight())
            return 0;
    }
    return 1;
}


int MySweepLDT::AddSharpTriangles(DoubleLinkedList* CurrP, int lrf)
// Method adds triangles left/right to the currently inserted triangle
// - CurrP: pointer in AF to the cirrently inersted triangle
// - lrf: Left Right Flag 
// returns pointer to the last nex position in AF 
// or NULL if the move has not been done
{
    DoubleLinkedList* AFWalker = nullptr;
    DoubleLinkedList* AFWalkerNext = nullptr;    

    if (SetWalkers(CurrP, AFWalker, AFWalkerNext, lrf) == 0)
        return 0;

    int EqualFlag = Equal(CurrP->key, AFWalker->key, RoundingTolerance);
        
    int TInd;
    double dotp;
 
    dotp = DotProduct(Vertices[CurrP->VertexNo].x, Vertices[CurrP->VertexNo].y,
            Vertices[AFWalker->VertexNo].x, Vertices[AFWalker->VertexNo].y,
            Vertices[AFWalkerNext->VertexNo].x, Vertices[AFWalkerNext->VertexNo].y);

//    if (lrf == RIGHT_WALK)
        if ((EqualFlag) && (dotp > 0))
            dotp = 1;


     if (dotp > 0)  // if angle is less than 90 degrees
     {
        TInd = CloseTriangle(CurrP, AFWalker, AFWalkerNext, lrf);
   
        CheckDT(NextTriangPosition, TInd);

        int tSecond = GetSecondTriangleIndex(TriangArray[NextTriangPosition], TInd);
        CheckDT(NextTriangPosition, tSecond);

        NextTriangPosition++;
        return 1;
    }
    return 0;
}




int MySweepLDT::CloseTriangle(DoubleLinkedList* p1, DoubleLinkedList* p2, DoubleLinkedList* p3, int lrf)
// Function close a triangle in the right direction of the 
// projection of NewVertexIndex
{
    int i = p1->VertexNo;
    int j = p2->VertexNo;
    int k = p3->VertexNo;
    int T1, T2;

    if (lrf == RIGHT_WALK)
    {
        T1 = p1->TriangleIndex;
        T2 = p2->TriangleIndex;

        InsertTriangle(i, k, j, NONE, T2, T1);

        SetNeighbour(T1, i, j, NextTriangPosition);
        SetNeighbour(T2, k, j, NextTriangPosition);

    }
    else // left walk
    {
        T1 = p2->TriangleIndex;
        T2 = p3->TriangleIndex;

        InsertTriangle(i, j, k, T1, T2, NONE);
        SetNeighbour(T1, i, j, NextTriangPosition);
        SetNeighbour(T2, k, j, NextTriangPosition);

    }
    SetAllAFFlags(NextTriangPosition);
    AF.DeleteElement(p2, NextTriangPosition);
    
    return T2;
}



void MySweepLDT::SetAllAFFlags(int TIndex)
{
    if (TriangArray[TIndex].ij != NONE)
        SetAFFlag(TriangArray[TIndex].ij);

    if (TriangArray[TIndex].jk != NONE)
        SetAFFlag(TriangArray[TIndex].jk);

    if (TriangArray[TIndex].ki != NONE)
        SetAFFlag(TriangArray[TIndex].ki);
}



int MySweepLDT::GetSecondTriangleIndex(BaseTriangle& t, int ind)
{
    int r = t.ij;
    if ((r != NONE) && (r != ind))
        return r;

    r = t.jk;
    if ((r != NONE) && (r != ind))
        return r;

    return t.ki;
}

/*

int MySweepLDT::SolveNearToVerticalAF(DoubleLinkedList* CurrP, int lrf)
{

    DoubleLinkedList* AFWalker = nullptr;
    DoubleLinkedList* AFWalkerNext = nullptr;

    if (SetWalkers(CurrP, AFWalker, AFWalkerNext, lrf) == 0)
        return 0;
    double dy = abs(CurrP->y - AFWalker->y);
    if (Equal(dy, 0, RoundingTolerance))
        return 0;

    double dx = abs(CurrP->key - AFWalker->key);
    double dxdy = dx / dy;
    
    if (dxdy >= 1) // the angle is greater than 45 degrees
        return 0;

    // otherwise, flatten the AF, if the angle is smaller than 180 degrees

    double sa = SignedArea(CurrP->key, CurrP->y, AFWalker->key, AFWalker->y, AFWalkerNext->key, AFWalkerNext->y, RoundingTolerance);

    if (Equal(sa, 0, RoundingTolerance)) return 0;

    if ((lrf == RIGHT_WALK) && (sa < 0)) return 0;
    if ((lrf == LEFT_WALK) && (sa > 0)) return 0;

    int TInd = CloseTriangle(CurrP, AFWalker, AFWalkerNext, lrf);
    

    CheckDT(NextTriangPosition, TInd);

    int tSecond = GetSecondTriangleIndex(TriangArray[NextTriangPosition], TInd);
    CheckDT(NextTriangPosition, tSecond);

    NextTriangPosition++;

    return 1;
}

*/


int MySweepLDT::SolveNearToVerticalAF(DoubleLinkedList* CurrP, int lrf)
{

    DoubleLinkedList* AFWalker = nullptr;
    DoubleLinkedList* AFWalkerNext = nullptr;

    if (SetWalkers(CurrP, AFWalker, AFWalkerNext, lrf) == 0)
        return 0;
    double dy = abs(CurrP->y - AFWalker->y);
    if (Equal(dy, 0, RoundingTolerance))
        return 0;

    double dx = abs(CurrP->key - AFWalker->key);
    double dxdy = dx / dy;

    if (dxdy >= 1) // the angle is greater than 45 degrees
        return 0;

    // otherwise, flatten the AF, if the angle is smaller than 180 degrees

    double sa = SignedArea(CurrP->key, CurrP->y, AFWalker->key, AFWalker->y, AFWalkerNext->key, AFWalkerNext->y, RoundingTolerance);

    if (Equal(sa, 0, RoundingTolerance)) 
        return Groove(AFWalker, AFWalkerNext, lrf);

    if ((lrf == RIGHT_WALK) && (sa < 0)) 
        return Groove(AFWalker, AFWalkerNext, lrf);
    if ((lrf == LEFT_WALK) && (sa > 0)) 
        return Groove(AFWalker, AFWalkerNext, lrf);

    int TInd = CloseTriangle(CurrP, AFWalker, AFWalkerNext, lrf);


    CheckDT(NextTriangPosition, TInd);

    int tSecond = GetSecondTriangleIndex(TriangArray[NextTriangPosition], TInd);
    CheckDT(NextTriangPosition, tSecond);

    NextTriangPosition++;

    return 1;
}

int  MySweepLDT::Groove(DoubleLinkedList* w1, DoubleLinkedList* w2, int lrf)
{
    DoubleLinkedList* w3;
    DoubleLinkedList* stopP;

    if (lrf == RIGHT_WALK)
    {
        w3 = w2->RPointer;
        stopP = AF.ReturnRight();
    }
    else
    {
        w3 = w2->LPointer;
        stopP = AF.ReturnLeft();
    }
    
    if (w3 == stopP)
        return 0;

    double sa = SignedArea(w1->key, w1->y, w2->key, w2->y, w3->key, w3->y, RoundingTolerance);

    if (Equal(sa, 0, RoundingTolerance))
        return 0;

    if ((lrf == RIGHT_WALK) && (sa < 0))
        return 0;
    if ((lrf == LEFT_WALK) && (sa > 0))
        return 0;

    int TInd = CloseTriangle(w1, w2, w3, lrf);

    CheckDT(NextTriangPosition, TInd);

    int tSecond = GetSecondTriangleIndex(TriangArray[NextTriangPosition], TInd);
    CheckDT(NextTriangPosition, tSecond);

    NextTriangPosition++;

    return 1;
}




void MySweepLDT::ReturnAllTriangles(std::vector<STLStyleTriangle>& T)
{
    STLStyleTriangle Tr;

    for (int u = 0; u < NextTriangPosition; u++)
    {

        Tr.x1 = Vertices[TriangArray[u].i].x;
        Tr.y1 = Vertices[TriangArray[u].i].y;
        Tr.x2 = Vertices[TriangArray[u].j].x;
        Tr.y2 = Vertices[TriangArray[u].j].y;
        Tr.x3 = Vertices[TriangArray[u].k].x;
        Tr.y3 = Vertices[TriangArray[u].k].y;

        T.push_back(Tr);
    }
}


int MySweepLDT::TIsNotArtfTriangle(BaseTriangle& bt)
{
    if ((bt.i == 0) || (bt.i == 1))  return 0;
    if ((bt.j == 0) || (bt.j == 1))  return 0;
    if ((bt.k == 0) || (bt.k == 1))  return 0;
    return 1;
}

void MySweepLDT::ReturnDTTriangles(std::vector<STLStyleTriangle>& T)
{
    STLStyleTriangle Tr;

    for (int u = 0; u < NextTriangPosition; u++)
    {
       if (TIsNotArtfTriangle(TriangArray[u]))
        {
            Tr.x1 = Vertices[TriangArray[u].i].x;
            Tr.y1 = Vertices[TriangArray[u].i].y;
            Tr.x2 = Vertices[TriangArray[u].j].x;
            Tr.y2 = Vertices[TriangArray[u].j].y;
            Tr.x3 = Vertices[TriangArray[u].k].x;
            Tr.y3 = Vertices[TriangArray[u].k].y;

            T.push_back(Tr);
        }
    }
}


// Finalizacija 
void MySweepLDT::Finalize()
{
    FillToConvexHull(true);
    ConstructBootomAF();
    FillToConvexHull(false);

    AF.ClearAF();
}

constexpr double RAD_TO_DEG = 180.0 / 3.14159265358979323846;

double _AngleInBetween(double x1, double y1, double x2, double y2, double x3, double y3) {
    double dx1 = x1 - x2, dy1 = y1 - y2;
    double dx2 = x3 - x2, dy2 = y3 - y2;
    double angle = atan2(dy1 * dx2 - dy2 * dx1,
                          dx1 * dx2 + dy1 * dy2) * RAD_TO_DEG;
    return angle < 0.0 ? angle + 360.0 : angle;
}

double AngleInBetween(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double dx = x1 - x2;
    double dy = y1 - y2;

    double angle_rad = atan2(dy, dx);

    double alpha1 = angle_rad * (180.0 / M_PI);

    if (alpha1 < 0)
        alpha1 += 360.0;
    
    dx = x3 - x2;
    dy = y3 - y2;

    angle_rad = atan2(dy, dx);

    double alpha2 = angle_rad * (180.0 / M_PI);

    if (alpha2 < 0)
        alpha2 += 360.0;

    return (alpha1 - alpha2);
}

int MySweepLDT::FillToConvexHull(bool UpperFlag)
{

    int ContinueFlag = 1;
    int TInd1, TInd2;
    DoubleLinkedList* w1;
    DoubleLinkedList* w2;
    DoubleLinkedList* w3;

    double x1, y1, x2, y2, x3, y3;
    double sa;

    int addedT = 0;

 //   AF.PrintAFCoordinates("d:\\AFCOORDStart.TXT");

    for (int q = 0; q < 2; q++)
    while (ContinueFlag)
    {
        ContinueFlag = 0;
        w1 = AF.ReturnLeft()->RPointer;
        if (w1 == nullptr) return 1;
        w2 = w1->RPointer;
        if (w2 == nullptr) return 1;
        w3 = w2->RPointer;
        
        if (w3 == nullptr)
            return 1;
        //         FILE* f;

        while (w3 != AF.ReturnRight())
        {
            x1 = Vertices[w1->VertexNo].x;
            y1 = Vertices[w1->VertexNo].y;

            x2 = Vertices[w2->VertexNo].x;
            y2 = Vertices[w2->VertexNo].y;

            x3 = Vertices[w3->VertexNo].x;
            y3 = Vertices[w3->VertexNo].y;

            
            if (UpperFlag)
                sa = SignedArea(x1, y1, x2, y2, x3, y3, RoundingTolerance);
            else
                sa = SignedArea(x3, y3, x2, y2, x1, y1, RoundingTolerance);

            if (UpperFlag)
            {
                double sa1 = AngleInBetween(x1, y1, x2, y2, x3, y3);
                if ((sa1 < 180) && (sa1 > 0))
                    sa = 1.0;
            }
   


            if (sa > 0)
            {

                if (!Equal(sa, 0, RoundingTolerance))
                {
                    TInd1 = w1->TriangleIndex;
                    TInd2 = w2->TriangleIndex;

                    // Orientation is clockwise
                    addedT++;
                    InsertTriangle(w1->VertexNo, w3->VertexNo, w2->VertexNo, NONE, TInd2, TInd1);
                    SetNeighbour(TInd1, w1->VertexNo, w2->VertexNo, NextTriangPosition);
                    SetNeighbour(TInd2, w2->VertexNo, w3->VertexNo, NextTriangPosition);

                    AF.DeleteElement(w2, TInd1);

                    w1->TriangleIndex = NextTriangPosition;

                    CheckDT(NextTriangPosition, TInd1);
                    CheckDT(NextTriangPosition, TInd2);
                    NextTriangPosition++;

                    ContinueFlag = 1;

                    w2 = w3;
                    w3 = w3->RPointer;
                }
                else
                {
                    w1 = w2;
                    w2 = w3;
                    w3 = w3->RPointer;
                }
            }
            else
            {
                w1 = w2;
                w2 = w3;
                w3 = w3->RPointer;
            }
        }
    }

    addedT++;

    //AF.PrintAFCoordinates("d:\\AFCOORD.TXT");
    
    return 1;
}



int MySweepLDT::TwoVerticesArtificial(int ti)
// method checks whether both vertices of triangle ti are artificial
{
    if ((TriangArray[ti].i < 2) && (TriangArray[ti].j < 2))
        return 1;

    if ((TriangArray[ti].j < 2) && (TriangArray[ti].k < 2))
        return 1;

    if ((TriangArray[ti].k < 2) && (TriangArray[ti].i < 2))
        return 1;

    return 0;
}



int  MySweepLDT::ReturnNeighbouringTriangle(int ti, int va, int vb)
{
    if (TriangArray[ti].i == va)
        if (TriangArray[ti].j == vb)
            return TriangArray[ti].ij;

    if (TriangArray[ti].i == vb)
        if (TriangArray[ti].j == va)
            return TriangArray[ti].ij;

    if (TriangArray[ti].j == va)
        if (TriangArray[ti].k == vb)
            return TriangArray[ti].jk;

    if (TriangArray[ti].j == vb)
        if (TriangArray[ti].k == va)
            return TriangArray[ti].jk;

    if (TriangArray[ti].k == va)
        if (TriangArray[ti].i == vb)
            return TriangArray[ti].ki;

    if (TriangArray[ti].k == vb)
        if (TriangArray[ti].i == va)
            return TriangArray[ti].ki;

    return NONE;
}

void MySweepLDT::ConstructBootomAF()
// v funkciji je treba postaviti še "border edge", da se ne bo izvajala legalizacija na umetnih trikotnikih in 
{
    AF.ClearAFInBetween();

    DoubleLinkedList* L = AF.ReturnLeft()->RPointer;
    DoubleLinkedList* R = AF.ReturnRight()->LPointer;

    if (L == R)
        L = AF.ReturnLeft();

    int Guard = L->VertexNo;

    int TIndex = R->TriangleIndex;
    TriangArray[TIndex].ArtificialFlag = true;
    int v1 = R->VertexNo;
    int v2 = 1;

    int JustTwoVertices = 0;

    if (TwoVerticesArtificial(TIndex))
    {
        TriangArray[TIndex].ArtificialFlag = true;

        v2 = 0;
        TIndex = ReturnNeighbouringTriangle(TIndex, v1, v2);
        if (TIndex != NONE)
            TriangArray[TIndex].ArtificialFlag = true;
        else
            JustTwoVertices = 1;
    }

    if (!JustTwoVertices)
    {
        int v3 = ThirdTriangleVertex(TIndex, v1, v2);
        int Tinner = ReturnNeighbouringTriangle(TIndex, v1, v3);
        TriangArray[TIndex].ArtificialFlag = true;
        SetNeighbour(Tinner, v1, v3, NONE);
        AF.InsertRightFrom(L, Vertices[v3].x, Vertices[v3].y, v3, Tinner);

        while (v3 != Guard)
        {
            v1 = v3;
            TIndex = ReturnNeighbouringTriangle(TIndex, v1, v2);
            TriangArray[TIndex].ArtificialFlag = true;

            if (TwoVerticesArtificial(TIndex))
            {
                v2 = 0;
                TIndex = ReturnNeighbouringTriangle(TIndex, v1, v2);
                TriangArray[TIndex].ArtificialFlag = true;
            }
            v3 = ThirdTriangleVertex(TIndex, v1, v2);
            Tinner = ReturnNeighbouringTriangle(TIndex, v1, v3);

            if (v3 != Guard)
                AF.InsertRightFrom(L, Vertices[v3].x, Vertices[v3].y, v3, Tinner);
            else
                L->TriangleIndex = Tinner;

            SetNeighbour(Tinner, v1, v3, NONE);
        }
        R->TriangleIndex = NONE;
    }
}



int MySweepLDT::ThirdTriangleVertex(int ti, int v1, int v2)
{
    if ((TriangArray[ti].i == v1) && (TriangArray[ti].j == v2))
        return TriangArray[ti].k;

    if ((TriangArray[ti].i == v2) && (TriangArray[ti].j == v1))
        return TriangArray[ti].k;

    if ((TriangArray[ti].j == v1) && (TriangArray[ti].k == v2))
        return TriangArray[ti].i;

    if ((TriangArray[ti].j == v2) && (TriangArray[ti].k == v1))
        return TriangArray[ti].i;

    if ((TriangArray[ti].k == v1) && (TriangArray[ti].i == v2))
        return TriangArray[ti].j;

    if ((TriangArray[ti].k == v2) && (TriangArray[ti].i == v1))
        return TriangArray[ti].j;

    return NONE;
}

int MySweepLDT::PrintTriangles(char* fname)
{
    FILE* f = fopen(fname, "wt");

    if (f == NULL)
        return 0;

    for (int i = 0; i < NextTriangPosition; i++)
    {
        fprintf(f, "%d: %d %d %d %d %d %d", i, TriangArray[i].i, TriangArray[i].j, TriangArray[i].k,
            TriangArray[i].ij, TriangArray[i].jk, TriangArray[i].ki);
        if (TriangArray[i].AFFlag == true)
            fprintf(f, " true\n");
        else
            fprintf(f, " false\n");
    }

    fclose(f);

    return 1;
}


int MySweepLDT::PrintNVertices(int n, char* fname)
{
    FILE* f = fopen(fname, "wt");

    if (f == NULL)
        return 0;

    if (n > Vertices.size())
        n = (int)Vertices.size();

    fprintf(f, "%d\n", n);
    for (int i = 2; i < n; i++)
        fprintf(f, "%lf, %lf\n", Vertices[i].x, Vertices[i].y);

    fclose(f);

    return 1;
}
