#include <iostream>
#include <set>
#include <vector>
#include "GeometryEntities/Vector.h"
#include "GeometryEntities/Matrix.h"
#include "GeometryEntities/Point.h"
#include "GeometryEntities/Line.h"
#include "GeometryEntities/Segment.h"
#include "GeometryEntities/Plane.h"

#include "Splines/BezierCurve.h"

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
  Vector<float, 3> c(std::vector<float>({1.f, 1.f, 1.f}));
  Vector<float, 3> d(p);
  Vector<float, 3> e(std::set<float>({1.f, 2.f, 3.f}));
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
  assert(FindRelationship(a, c) == VectorRelationship::None);
  assert(FindRelationship(a, b) == VectorRelationship::Orthogonal);
  assert(a != b);
  assert(a == a);
}

void MatrixTest() {
  /// constructors
  Matrix<double, 3> a({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
  Matrix<double, 3> b(Vector3d(1, 0, 1), Vector3d(0, 1, 2), Vector3d(3, 1, 0));
  Matrix<double, 3, 1> c({{1, 1, 1}});
  Matrix<double, 3> e;
  /// arithmetic
  std::cout << a << '\n' << b << '\n' << e << '\n';
  std::cout << a + b << '\n' << a - b << '\n';
  std::cout << a * b << '\n';
  std::cout << a * c << '\n';
  auto copy = a;
  copy *= b;
  std::cout << copy << '\n';
  std::cout << a + 3 << '\n' << a - 5 << '\n' << a * 6 << '\n';

  /// matrix calc
  std::cout << Transposed(c) << '\n';
  std::cout << Triangulated(a) << '\n';
  std::cout << Diagonalized(b) << '\n';
  assert(a.Determinant() == 0);
  assert(b.Determinant() == -5);
  assert(a.Trace() == 15);
  assert(b.Trace() == 2);
  assert(a.Rank() == 2);
  assert(b.Rank() == 3);
  assert(c.Rank() == 1);
  assert(b * Inverted(b) == e);
  /// relationship
  assert(a == a);
  assert(a != b);
}

void PointTest() {
  /// constructors
  Point3f a(1.f, 0.f, 0.f);
  Point3<float> b({0.f, 1.f, 0.f});
  Point<float, 3> c(vector<float>({1.f, 1.f, 1.f}));
  cout << a << " " << b << " " << c << '\n';
  /// arithmetic
  cout << a + b << " " << a - b << " " << c * 5 << " " << 3 * c / 5 << " " << a * c << '\n';

  /// calc
  cout << DotProduct(a, c) << " " << DotProduct(a, b) << " " << CrossProduct(a, b) << " " << MixedProduct(a, b, c)
       << " " << a.SquaredDistance(b) << " " << a.Distance(c) << " " << Normalised(c) << " " << Cos(a, b) << " "
       << Sin(a, c) << " " << Angle(a, b) << '\n';

  /// relationship
  assert(FindRelationship(a, c) == PointRelationship::None);
  assert(FindRelationship(a, a) == PointRelationship::Identical);
  assert(a != b);
  assert(a == a);
}

void LineTest() {
  /// constructors
  Line3d a(Point3d(0, 0, 0), Point3d(1, 1, 1));
  Line3d b(Point3d(0, 0, 0), Vector3d(1, 0, 0));
  Line3d c;
  c.GetOrigin() = {3, 1, 2};
  c.GetDirection() = {1, 1, 1};
  Line3d d(a.GetPoint(10), Normalised(a.GetDirection()));
  std::cout << "constructors ok\n";
  /// setters and getters
  assert(a[10] == Point3d(10, 10, 10));
  assert(c[-1] == Point3d(2, 0, 1));
  std::cout << "getters ok\n";
  /// calc
  auto copy = a;
  copy.Normalise();
  assert(copy.GetDirection().Length() == 1.);
  assert(copy.GetPoint(10 * a.GetDirection().Length()) == a.GetPoint(10));
  assert(b.Contains({5, 0, 0}));
  assert(!b.Contains({5, 5, 5}));
  std::cout << "Contains ok\n";

  assert(a.SquaredDistance({0, 0, 0}) == 0);
  assert(b.SquaredDistance({0, 0, 1}) == 1);

  auto copyb = b;
  copyb.Normalise();
  copyb.GetOrigin() = {0, 0, 1};
  assert(b.SquaredDistance(copyb) == 1);
  assert(b.SquaredDistance(a) == 0);
  assert(b.SquaredDistance(b) == 0);
  std::cout << "distance ok\n";

  auto inter = a.Intersection(b);
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == a.GetOrigin());

  inter = a.Intersection(a);
  assert(inter->GetType() == Entity::Line);
  assert(*static_cast<Line3d*>(inter.get()) == a);

  inter = a.Intersection(Segment3d(a.GetPoint(0), a.GetPoint(1)));
  assert(inter->GetType() == Entity::Segment);
  assert(*static_cast<Segment3d*>(inter.get()) == Segment3d(a.GetPoint(0), a.GetPoint(1)));

  inter = a.Intersection(Segment3d(b.GetPoint(-1), b.GetPoint(1)));
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == a.GetOrigin());

  inter = a.Intersection(Segment3d(b.GetPoint(-1), b.GetPoint(-0.5)));
  assert(inter->GetType() == Entity::Void);

  std::cout << "intersection ok\n";

  assert(b.Projection({0, 1, 0}) == Point3d(0, 0, 0));
  assert(a.Projection({1, 1, 1}) == Point3d(1, 1, 1));

  assert(b.Projection(Segment3d({-1, 0, 0}, {1, 0, 0})) == Segment3d({-1, 0, 0}, {1, 0, 0}));
  assert(b.Projection(Segment3d({-1, 0, 5}, {1, 0, 10})) == Segment3d({-1, 0, 0}, {1, 0, 0}));

  std::cout << "projection ok\n";
  /// relationship
  assert(FindRelationship(a, c) == LineRelationship::Parallel);
  assert(FindRelationship(a, a) == LineRelationship::Identical);
  assert(FindRelationship(a, b) == LineRelationship::Intersecting);
  assert(FindRelationship(a, copyb) == LineRelationship::Skew);
  assert(a != b);
  assert(a == d);
}

void SegmentTest() {
  /// constructors
  Segment3d a(Point3d(-3, 0, 0), Point3d(1, 0, 0));
  Segment3d b(Point3d(0, 0, 0), Point3d(1, 0, 0));
  Segment3d c(Point3d(0, -1, 0), Vector3d(0, 1, 0));
  Segment3d d(Point3d(2, 0, 0), Vector3d(2, 1, 0));
  std::cout << "constructors ok\n";
  /// setters and getters
  assert(a.GetPoint(0.75) == Point3d(0, 0, 0));
  assert(a[1] == a.GetRight());
  assert(a.GetLeft() == b.GetPoint(-3));
  assert(a.GetDirection() == Vector3d(4, 0, 0));
  assert(a.GetMidpoint() == Point3d(-1, 0, 0));
  std::cout << "getters ok\n";

  /// calc
  assert(a.Contains(Point3d(0, 0, 0)));
  assert(!a.Contains(Point3d(1, 1, 1)));

  std::cout << "Contains ok\n";

  assert(a.SquaredDistance(Point3d(0, 1, 0)) == 1);
  assert(a.SquaredDistance(Point3d(0, 0, 0)) == 0);
  assert(a.SquaredDistance(Point3d(2, 1, 0)) == 2);

  assert(a.SquaredDistance(b) == 0);
  assert(a.SquaredDistance(c) == 0);
  assert(c.SquaredDistance(d) == 4);
  assert(a.SquaredDistance(d) == 1);
  assert(a.SquaredDistance(Segment3d({2, 0, 0}, {3, 1, 1})) == 1);

  std::cout << "distance ok\n";

  auto inter = a.Intersection(Point3d(0, 0, 0));
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == Point3d(0, 0, 0));

  inter = a.Intersection(b);
  assert(inter->GetType() == Entity::Segment);
  assert(*static_cast<Segment3d*>(inter.get()) == b);

  inter = a.Intersection(c);
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == Point3d(0, 0, 0));

  inter = a.Intersection(d);
  assert(inter->GetType() == Entity::Void);
  std::cout << "intersection ok\n";

  assert(a.Projection(Point3d(0, 15, 2)) == Point3d(0, 0, 0));
  assert(c.Projection(d) == Segment3d({0, 0, 0}, {0, 1, 0}));

  std::cout << "projection ok\n";

  /// relationship
  assert(FindRelationship(a, c) == SegmentRelationship::Intersecting);
  assert(FindRelationship(a, b) == SegmentRelationship::Intersecting);
  assert(FindRelationship(a, a) == SegmentRelationship::Identical);
  assert(FindRelationship(a, d) == SegmentRelationship::Skew);
  assert(FindRelationship(d, c) == SegmentRelationship::Parallel);
  assert(a != b);
  assert(a == a);
}

void PlaneTest() {
  /// constructors
  Plane<double> a(Point3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(0, 2, 0));
  Plane<double> b(Point3d(1, 0, 1), Point3d(2, 0, 1), Point3d(0, 1, 1));
  Plane<double> c(Point3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(0, 0, 1));
  Plane<double> d;
  std::cout << "constructors ok\n";

  /// setters and getters
  d.GetOrigin() = {0, 0, 0};
  d.GetAbscissa() = {0, 1, 0};
  d.GetOrdinate() = {0, 0, 1};
  assert(a.GetPoint(2, 1) == Point3d(2, 2, 0));
  assert(a.GetNormal() == Vector3d(0, 0, 2));
  std::cout << "getters ok\n";

  /// calc
  auto copy = a;
  copy.Normalise();
  assert(a.Contains(Point3d(1, 1, 0)));
  assert(!a.Contains(Point3d(1, 1, 10)));
  assert(a.Contains(Line3d(Point3d{0, 1, 0}, Point3d{1, 0, 0})));
  assert(a.Contains(Segment3d({1, 1, 0}, {0, 1, 0})));
  assert(!a.Contains(Segment3d({1, 1, 0}, {0, 1, 1})));

  std::cout << "Contains ok\n";

  assert(a.SquaredDistance(b) == 1);
  assert(a.SquaredDistance(Point3d(1, 15, -3)) == 9);
  assert(a.SquaredDistance(Line3d({0, 0, 0}, Point3d(1, 1, 1))) == 0);
  assert(a.SquaredDistance(Line3d({1, 1, 1}, Point3d(1, 0, 1))) == 1);
  assert(a.SquaredDistance(Segment3d({100, 1, 1}, Point3d(1, 0, 2))) == 1);

  std::cout << "distance ok\n";

  auto inter = a.Intersection(Point3d(1, 1, 0));
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == Point3d(1, 1, 0));

  inter = a.Intersection(Line3d({0, 0, -1}, Point3d(0, 0, 1)));
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == Point3d(0, 0, 0));

  inter = a.Intersection(Line3d({0, 1, 0}, Point3d(1, 0, 0)));
  assert(inter->GetType() == Entity::Line);
  assert(*static_cast<Line3d*>(inter.get()) == Line3d({0, 1, 0}, Point3d(1, 0, 0)));

  inter = a.Intersection(Segment3d({0, 1, 1}, Point3d(1, 0, 0)));
  assert(inter->GetType() == Entity::Point);
  assert(*static_cast<Point3d*>(inter.get()) == Point3d(1, 0, 0));

  inter = a.Intersection(Segment3d({0, 1, 0}, Point3d(1, 0, 0)));
  assert(inter->GetType() == Entity::Segment);
  assert(*static_cast<Segment3d*>(inter.get()) == Segment3d({0, 1, 0}, Point3d(1, 0, 0)));

  inter = a.Intersection(c);
  assert(inter->GetType() == Entity::Line);
  assert(*static_cast<Line3d*>(inter.get()) == Line3d(Point3d(0, 0, 0), Vector3d(1, 0, 0)));

  inter = a.Intersection(copy);
  assert(inter->GetType() == Entity::Plane);
  assert(*static_cast<Plane<double>*>(inter.get()) == a);
  std::cout << "intersection ok\n";

  assert(a.Projection(Point3d(1, 1, 1)) == Point3d(1, 1, 0));
  assert(a.Projection(Line3d({1, 0, 1}, Point3d(0, 1, 0))) == Line3d({1, 0, 0}, Point3d(0, 1, 0)));
  assert(a.Projection(Segment3d({1, 2, 3}, {4, 5, 6})) == Segment3d({1, 2, 0}, {4, 5, 0}));

  std::cout << "projection ok\n";

  /// relationship
  assert(FindRelationship(a, c) == PlaneRelationship::Intersecting);
  assert(FindRelationship(a, b) == PlaneRelationship::Parallel);
  assert(FindRelationship(a, copy) == PlaneRelationship::Identical);
  assert(a != b);
  assert(a == a);
}

void GeometryTest() {
  cout << "Vector Test \n";
  VectorTest();

  cout << "\nMatrix Test \n";
  MatrixTest();

  cout << "\nPoint Test \n";
  PointTest();

  cout << "\nLine Test \n";
  LineTest();

  cout << "\nSegment Test \n";
  SegmentTest();

  cout << "\nPlane Test \n";
  PlaneTest();
}

void BezierTest() {

}

int main() {
  cout << "\nGeometry Test \n";
  GeometryTest();
}
