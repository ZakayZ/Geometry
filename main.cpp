#include <iostream>
#include "Geometry/Vector.h"
#include "Geometry/Point.h"
#include "Geometry/Line.h"
#include "Geometry/Segment.h"
#include "Geometry/Plane.h"

using namespace std;

void VectorTest() {
  /// Data
  Data<float, 3> p;
  p.x = 1;
  p.y = 2;
  p.z = 3;


  /// constructors
  Vector<float, 3> a(1.f, 0.f, 0.f);
  Vector<float, 3> b({0.f, 1.f, 0.f});
  Vector<float, 3> c(vector<float>({1.f, 1.f, 1.f}));
  Vector<float, 3> d(p);
  cout << a << " " << b << " " << c << " " << d << '\n';
  d += d;
  d *= 5;
  d /= 5;
  d -= d;
  /// arithmetic
  cout << a + b << " " << a - b << " " << c * 5 << " " << 3. * c / 5. << " " << a * c << '\n';

  /// vector calc
  cout << DotProduct(a, c) << " " << DotProduct(a, b) << " " << CrossProduct(a, b) << " " << MixedProduct(a, b, c)
       << " " << a.SquaredLength() << " " << b.Length() << " " << Normalised(c) << " " << Cos(a, b) << " "
       << Sin(a, c) << " " << Angle(a, b) << '\n';

  /// relationship
  assert(FindRelationShip(a, c) == Vector3f::Relationship::None);
  assert(FindRelationShip(a, b) == Vector3f::Relationship::Orthogonal);
  assert(a != b);
  assert(a == a);
}

void PointTest() {
  /// constructors
  Point3f a(1.f, 0.f, 0.f);
  Point3<float> b({0.f, 1.f, 0.f});
  Point<float, 3> c(vector<float>({1.f, 1.f, 1.f}));
  cout << a << " " << b << " " << c << '\n';
  /// arithmetic
  cout << a + b << " " << a - b << " " << c * 5 << " " << 3 * c / 5 << " " << a * c << '\n';

  /// vector calc
  cout << DotProduct(a, c) << " " << DotProduct(a, b) << " " << CrossProduct(a, b) << " " << MixedProduct(a, b, c)
       << " " << SquaredDistance(a, b) << " " << Distance(a, c) << " " << Normalised(c) << " " << Cos(a, b) << " "
       << Sin(a, c) << " " << Angle(a, b);

  /// relationship
  assert(FindRelationShip(a, c) == Point3f::Relationship::None);
  assert(FindRelationShip(a, a) == Point3f::Relationship::Identical);
  assert(a != b);
  assert(a == a);
}

int main() {
  cout << "Vector Test \n";
  VectorTest();

  cout << "\nPoint Test \n";
  PointTest();
}
