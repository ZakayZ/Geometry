#include <iostream>
#include "Geometry/Vector.h"
#include "Geometry/Line.h"

using namespace std;

int main() {
  Point3d p;
  p.x = 1;
  p.y = 2;
  p.z = 3;
  Vector<float, 3> a(1.f, 0.f, 0.f);
  Vector<float, 3> b(0.f, 1.f, 0.f);
  cout << SquaredDistance(a, b);

  Line<double, 3> line(a, b);
}
