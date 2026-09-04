
#include "AuxF.h"



#include <limits>
#include <cmath>
#include <math.h>


double MaxDouble() { return (std::numeric_limits<double>::max)(); }



int Equal(double a, double b, double delta)
{
	if (std::abs(a - b) <= delta) return 1; else return 0;
}


int EqualPoints(double x1, double y1, double x2, double y2, double delta)
{
	if (Equal(x1, x2, delta) && Equal(y1, y2, delta))
		return 1;
	return 0;
}



double Pitagora(double x1, double y1, double x2, double y2)
{
	double a = (x2 - x1) * (x2 - x1);
	double b = (y2 - y1) * (y2 - y1);
	return sqrt(a + b);
}

double Perimeter(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a = Pitagora(x1, y1, x2, y2);
	double b = Pitagora(x2, y2, x3, y3);
	double c = Pitagora(x3, y3, x1, y1);
	return (a + b + c);

}


int FromRealToPixelX(double x, double xwmin, double dx, int pixelsize)
// converts real coordinate to pixel coordinate and returns it
// - x: coordinate which will be converted
// - xwmin: minimal value of the real window
// - dx: size of a real window
// - pixelsize: size of the window in pixels
{
	return ((int)((double)((x - xwmin) * pixelsize) / dx));
}

int FromRealToPixelY(double y, double ywmin, double dy, int pixelsize)
// converts real coordinate to pixel coordinate and returns it
// - y: coordinate which will be converted
// - ywmin: minimal value of the real window
// - dy: size of a real window
// - pixelsize: size of the window in pixels
{
	return (FromRealToPixelX(y, ywmin, dy, pixelsize));

	//  return ( pixelsize - ((int) ( (double)((y - ywmin) * pixelsize) / dy)) );
}


int FromRealToPixel(double x, double xwmin, double dx, int pixelsize)
// converts real coordinate to pixel coordinate and returns it
// - x: coordinate which will be converted
// - xwmin: minimal value of the real window
// - dx: size of a real window
// - pixelsize: size of the window in pixels
// The function can be aplied for X and Y direction
{
	return ((int)((float)((x - xwmin) * pixelsize) / dx));
}



int TheSamePoints2D(double x1, double y1, double x2, double y2, double delta)
{
	if ((Equal(x1, x2, delta) == 1) && (Equal(y1, y2, delta) == 1)) return 1;
	return 0;
}


int TheSamePoints2D(double x1, double y1, double x2, double y2)
{
	return TheSamePoints2D(x1, y1, x2, y2, MYEPSILON);
}


void Normalise(double& x, double& y)
{
	double length = sqrt(x * x + y * y);
	if (length == 0)
	{
		x = y = 0.0;
	}
	else
	{
		x /= length;
		y /= length;
	}
}


double DotProduct(double x1, double y1, double x2, double y2, double x3, double y3)
{
//	Normalise(x1, y1);
//	Normalise(x2, y2);
//	Normalise(x3, y3);

	double ax = x1 - x2;
	double bx = x3 - x2;
	double ay = y1 - y2;
	double by = y3 - y2;

	return(ax * bx + ay * by);
}


double SignedArea(double x1, double y1, double x2, double y2, double x3, double y3, double delta)
{
	double a;

//	Normalise(x1, y1);
//	Normalise(x2, y2);
//	Normalise(x3, y3);


	if (Equal(y1, y2, delta) && Equal(y2, y3, delta))
		return 0;



	a = x1 * (y2 - y3) - y1 * (x2 - x3) + (x2 * y3 - y2 * x3);
	return a;
}


double SignedAreaNormalized(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a;

	Normalise(x1, y1);
	Normalise(x2, y2);
	Normalise(x3, y3);

	a = x1 * (y2 - y3) - y1 * (x2 - x3) + (x2 * y3 - y2 * x3);
	return a;
}