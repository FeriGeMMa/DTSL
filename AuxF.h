#pragma once


#define PI 3.141592654
#define PI2 6.283185307
#define PI_HALF 1.570796327
#define MYEPSILON 0.0000001
#define GAP_BORDER 1.0

#define NONE -1


double MaxDouble();


int Equal(double a, double b, double delta);


int FromRealToPixelX(double x, double xwmin, double dx, int pixelsize);
int FromRealToPixelY(double y, double ywmin, double dy, int pixelsize);
int FromRealToPixel(double x, double xwmin, double dx, int pixelsize);


int TheSamePoints2D(double x1, double y1, double x2, double y2, double delta);
int TheSamePoints2D(double x1, double y1, double x2, double y2);


void Normalise(double& x, double& y);


double DotProduct(double x1, double y1, double x2, double y2, double x3, double y3);
double SignedArea(double x1, double y1, double x2, double y2, double x3, double y3, double delta);
double SignedAreaNormalized(double x1, double y1, double x2, double y2, double x3, double y3);

double Pitagora(double x1, double y1, double x2, double y2);
double Perimeter(double x1, double y1, double x2, double y2, double x3, double y3);

// just for texting the performance
extern int AFSearchCounter;
extern int AFJumpCounter;