#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "MySweepLDT.h"

int exampleUsage() {
  std::ifstream file("samples/clus1000000.pnt");

  if (!file) {
    std::cerr << "Cannot open points.txt\n";
    return 1;
  }

  double BoxXMax, BoxYMax, BoxXMin, BoxYMin;
  std::vector<double> XArr, YArr;

  BoxXMin = BoxYMin = MaxDouble();
  BoxXMax = BoxYMax = -MaxDouble();

  std::string line;

  int NoOfPoints = 0;
  bool f = true;
  int i = 0;
  while (std::getline(file, line)) {
    if (line.empty()) continue;

    std::stringstream ss(line);

    if (f) {
      f = false;

      ss >> NoOfPoints;

      XArr.reserve(NoOfPoints);
      YArr.reserve(NoOfPoints);

      continue;
    }

    double x, y;
    char comma;

    if (ss >> x >> comma >> y && comma == ',') {
      XArr.push_back(x);
      YArr.push_back(y);

      if (XArr[i] > BoxXMax) BoxXMax = XArr[i];
      if (YArr[i] > BoxYMax) BoxYMax = YArr[i];
      if (XArr[i] < BoxXMin) BoxXMin = XArr[i];
      if (YArr[i] < BoxYMin) BoxYMin = YArr[i];

      i++;
    }
  }

  auto start = std::chrono::high_resolution_clock::now();

  MySweepLDT DT;
  std::vector<STLStyleTriangle> DTTriangles;

  int LTEntries = (int)XArr.size();
  LTEntries = (int)(LTEntries / 100);
  if (LTEntries == 0) LTEntries = 2;

  int DataFlag =
      DT.InitSL(BoxXMin, BoxYMin, BoxXMax, BoxYMax, LTEntries, XArr, YArr);

  if (DataFlag == 0) return -1;

  DT.DoDT();

  auto end = std::chrono::high_resolution_clock::now();

  double milliseconds =
      std::chrono::duration<double, std::milli>(end - start).count();

  std::cout << "Triangulation time: " << milliseconds << " ms\n";
  std::cout << "Triangles: " << DTTriangles.size() << '\n';

  return -1;
}

int main() { exampleUsage(); }
