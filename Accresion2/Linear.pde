/*
SIMPLE RAY TRACING SKETCH
 */
double colorTable[];
final double MAX_MAP=9e16f;
final double E=0.00001f;
final double Epsilon = 0.001f;
final double PI90 = (double)Math.PI * 0.5f;
final double PI180 = (double)Math.PI;
final double PI2 = (double)Math.PI * 2;
final double InvPI2 = (double)(0.5 / Math.PI);
final double InvPI = (double)(1.0 / Math.PI);
final double PI270 = PI90 * 3.0f;
final double PISqrd=9.8696044010893586188344909998762;
    
final byte PositionBetween=0, PositionBeyond=1, PositionOrigin=2, PositionDestination=3, PositionBehind=4, PositionAbove=5, PositionBelow=6, PositionLeft=7, PositionRight=8;
final byte AxisX=0, AxisY=1, AxisZ=2, AxisW=3; 
final byte SideLeft=1, SideRight=2, SideTop=4, SideBottom=8, SideFront=16, SideBack=32;
class Linear
{
  public int x;
  public int y;
  public int z;
}
class Local extends Linear
{
  Local() { 
    x = 0; 
    y = 0; 
    z = 0;
  }
  Local(short xx, short yy) { 
    x = xx; 
    y = yy; 
    z = 0;
  }
  Local(int xx, int yy) { 
    x = xx; 
    y = yy; 
    z = 0;
  }
  Local(double xx, double yy) { 
    x = (int)xx; 
    y = (int)yy; 
    z = 0;
  }
}       
boolean Between2D(Cartesian minimum, Cartesian at, Cartesian maximum)
{
  return (minimum.x <= at.x && at.x <= maximum.x) && (minimum.y <= at.y && at.y <= maximum.y);
}
boolean Between3D(Cartesian minimum, Cartesian at, Cartesian maximum)
{
  return (minimum.x <= at.x && at.x <= maximum.x) && (minimum.y <= at.y && at.y <= maximum.y) && (minimum.z <= at.z && at.z <= maximum.z);
}
boolean Overlapping2D(Cartesian minimum, Cartesian maximum)
{
  return (minimum.x <= maximum.x) && (minimum.y <= maximum.y);
}
boolean Overlapping2D(Cartesian minimum, Cartesian maximum, double threshold)
{
  return (minimum.x <= maximum.x + threshold) && (minimum.y <= maximum.y + threshold);
}
boolean Overlapping3D(Cartesian minimum, Cartesian maximum)
{
  return (minimum.x <= maximum.x) && (minimum.y <= maximum.y) && (minimum.z <= maximum.z);
}
boolean Overlapping3D(Cartesian minimum, Cartesian maximum, double threshold)
{
  return (minimum.x <= maximum.x + threshold) && (minimum.y <= maximum.y + threshold) && (minimum.z <= maximum.z + threshold);
}
double Clamp(double v, double min, double max)
{
  if (v < min) { 
    v = min;
  } else if (v > max) { 
    v = max;
  }
  return v;
}
void Clamp(Cartesian c){
  c.x=Clamp(c.x,0,1);
  c.y=Clamp(c.y,0,1);
  c.z=Clamp(c.z,0,1);
}
class Point extends Cartesian {

  double r;
  Point(double x, double y) {
    super(x, y);
    r=1;
  }
  Point(double rad) {
    r=rad;
  }
  Point() {
    r=1;
  }
  Point(double x, double y, double z) {
    super(x, y, z);
    r=1;
  }
  Point(double x, double y, double z, double radius) {
    super(x, y, z);
    r=radius;
  }
  Vector Sub(Point p1, Point p2) {
    return new Vector(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z);
  }
  boolean Parallel3D(Point p, double r)
  {
    Vector na = CrossProduct(this, p).AsVector();
    return r * r < Math.abs(DotProduct(na, p));
  }
  double Sine(double invr)
  {
    return x * invr;
  }
  double Cosine(double invr)
  {
    return y * invr;
  }
  Vector Absolute()
  {
    return new Vector(Math.abs(x), Math.abs(y), Math.abs(z));
  }
  Point Scaled2D(double s)
  {

    return new Point(x * s, y * s);
  }
  Point Scaled3D(double s)
  {

    return new Point(x * s, y * s, z * s);
  }
  Point Scaled3D(Vector s)
  {

    return new Point(x * s.x, y * s.y, z * s.z);
  }
  void Scale3D(Vector s)
  {

    x *= s.x;
    y *= s.y;
    z *= s.z;
  }
  void Scale2D(Vector s)
  {

    x *= s.x;
    y *= s.y;
  }
  Local Localize(double s)
  {
    short xs = (short)(x * s);
    short ys = (short)(y * s);
    return new Local(xs, ys);
  }
  Point Localize(Local loc, double s)
  {
    double iss = 1.0f / s;
    double xs = (double)(loc.x) * iss;
    double ys = (double)(loc.y) * iss;
    return new Point(xs, ys);
  }
  Point Centralize(Local loc, double s)
  {
    double iss = 1.0f / s;
    double xs = (((double)(loc.x) + 0.5f) * iss);
    double ys = (((double)(loc.y) + 0.5f) * iss);
    return new Point(xs, ys);
  }
  Point Localize(Local loc, double h, double s)
  {
    double iss = 1.0f / s;
    double xs = (double)(loc.x) * iss;
    double ys = (double)(loc.y) * iss;
    return new Point(xs, h, ys);
  }
  Point Centralize(Local loc, double h, double s)
  {
    double iss = 1.0f / s;
    double xs = (((double)(loc.x) + 0.5f) * iss);
    double ys = (((double)(loc.y) + 0.5f) * iss);
    return new Point(xs, h, ys);
  }
  double Max()
  {
    double mxy = Math.max(x, y);
    double mxz = Math.max(mxy, z);
    return mxz;
  }
  double Min()
  {
    double mxy = Math.min(x, y);
    double mxz = Math.min(mxy, z);
    return mxz;
  }

  Point Max(Point m)
  {
    return new Point(Math.max(x, m.x), Math.max(y, m.y), Math.max(z, m.z));
  }
  Point Min(Point m)
  {
    return new Point(Math.min(x, m.x), Math.min(y, m.y), Math.min(z, m.z));
  }

  Point Set3D(double xyz[])
  {
    x = xyz[0];
    y = xyz[1];
    z = xyz[2];
    return this;
  }
  Point Set2D(double xyz[])
  {
    x = xyz[0];
    y = xyz[1];
    return this;
  }
  boolean EqualSphere(Point p, double r)
  {
    Vector a = Sub(this, p);
    double r2 = (this.r+r) *(this.r* r);
    if (a.LengthSquared() < r2) return true;
    return false;
  }
  boolean EqualCube(Point p, double threshold)
  {
    return Math.abs(x - p.x) < threshold && Math.abs(y - p.y) < threshold && Math.abs(z - p.z) < threshold;
  }
  double Distance( Point p2)
  {
    return Sub(this, p2).Length();
  }
  double DistanceSquared(Point p2)
  {
    return Sub(this, p2).LengthSquared();
  }
  int PointTest(Point p, Point q)
  {
    if (LessThan(p, q))
    {
      return -1;
    }
    if (Greater(p, q))
    {
      return 1;
    }
    return 0;
  }
  int PointTest( Point q)
  {
    if (LessThan(this, q))
    {
      return -1;
    }
    if (Greater(this, q))
    {
      return 1;
    }
    return 0;
  }
  int TestPoints(Point p, Point q)
  {
    if (LessThan(p, q))
    {
      return -1;
    }
    if (Greater(p, q))
    {
      return 1;
    }
    return 0;
  }
  int PointAreaTest(Point p, Point q)
  {
    double pointAreaWidth=p.r+q.r;
    double np = p.x + p.y * pointAreaWidth;
    double nq = q.x + q.y * pointAreaWidth;
    if (np < nq)
    {
      return -1;
    }
    if (np > nq)
    {
      return 1;
    }
    return 0;
  }
  int SideCompare2D(Point a, Point b, byte side)
  {
    switch (side)
    {
    case SideLeft:
      {

        if (a.x < b.x)
        {
          return -1;
        }
        if (a.x > b.x)
        {
          return 1;
        }
        break;
      }
    case SideRight:
      {

        if (a.x > b.x)
        {
          return -1;
        }
        if (a.x < b.x)
        {
          return 1;
        }
        break;
      }
    case SideFront:
      {

        if (a.y < b.y)
        {
          return -1;
        }
        if (a.y > b.y)
        {
          return 1;
        }
        break;
      }
    case SideBack:
      {

        if (a.y > b.y)
        {
          return -1;
        }
        if (a.y < b.y)
        {
          return 1;
        }
        break;
      }
    case SideTop:
      {

        if (a.z < b.z)
        {
          return -1;
        }
        if (a.z > b.z)
        {
          return 1;
        }
        break;
      }
    case SideBottom:
      {

        if (a.z > b.z)
        {
          return -1;
        }
        if (a.z < b.z)
        {
          return 1;
        }
        break;
      }
    }
    return 0;
  }
  int SideCompare3D(Point a, Point b, byte side)
  {
    switch (side)
    {
    case SideLeft:
      {

        if (a.x < b.x)
        {
          return -1;
        }
        if (a.x > b.x)
        {
          return 1;
        }
        break;
      }
    case SideRight:
      {

        if (a.x > b.x)
        {
          return -1;
        }
        if (a.x < b.x)
        {
          return 1;
        }
        break;
      }
    case SideFront:
      {

        if (a.z < b.z)
        {
          return -1;
        }
        if (a.z > b.z)
        {
          return 1;
        }
        break;
      }
    case SideBack:
      {

        if (a.z > b.z)
        {
          return -1;
        }
        if (a.z < b.z)
        {
          return 1;
        }
        break;
      }
    case SideTop:
      {

        if (a.y < b.y)
        {
          return -1;
        }
        if (a.y > b.y)
        {
          return 1;
        }
        break;
      }
    case SideBottom:
      {

        if (a.y > b.y)
        {
          return -1;
        }
        if (a.y < b.y)
        {
          return 1;
        }
        break;
      }
    }
    return 0;
  }
  int LeftToRightCompare(Point a, Point b)
  {
    if (LessThan(a, b))
    {
      return -1;
    }
    if (Greater(a, b))
    {
      return 1;
    }
    return 0;
  }
  int RightToLeftCompare(Point a, Point b)
  {
    if (LessThan(b, a))
    {
      return -1;
    }
    if (Greater(b, a))
    {
      return 1;
    }
    return 0;
  }
}
class  Normal extends Cartesian {
  Normal(double x, double y, double z) {
    super(x, y, z);
  }
  Normal(double x, double y) {
    super(x, y);
  }
}
class Cartesian {

  //
  double x;
  double y;
  double z;
  Cartesian()
  {
    this.x = 0;
    this.y = 0;
    this.z = 0;
  }
  Cartesian(double x, double y)
  {
    this.x = x;
    this.y = y;
    this.z = 0;
  }
  Cartesian(double x, double y, double z)
  {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  Cartesian(Point p)
  {
    this.x = p.x;
    this.y = p.y;
    this.z = p.z;
  }
  Cartesian(Vector p)
  {
    this.x = p.x;
    this.y = p.y;
    this.z = p.z;
  }
  Cartesian(Normal p)
  {
    this.x = p.x;
    this.y = p.y;
    this.z = p.z;
  }
  int Orientation2D(Point p0, Point p1, Point p2)
  {
    Vector a = p1.Sub(p1, p0);
    Vector b = p2.Sub(p2, p0);
    double sa = a.x * b.y - a.y * b.x;
    if (sa > 0)
    {
      return 1;
    }
    if (sa < 0)
    {

      return -1;
    }
    return 0;
  }
  int Orientation2D(Cartesian p1, Cartesian p2)
  {
    Cartesian a = Sub(p1, this);
    Cartesian b = Sub(p2, this);
    double sa = a.x * b.y - a.y * b.x;
    if (sa > 0)
    {
      return 1;
    }
    if (sa < 0)
    {

      return -1;
    }
    return 0;
  }
  int Orientation3D(Cartesian p0, Cartesian p1, Cartesian p2)
  {
    Cartesian a = Sub(p1, p0);
    Cartesian b = Sub(p2, p0);
    Cartesian sa =a.CrossProduct(a, b);
    double ln = sa.x + sa.y + sa.z;
    if (ln > 0)
    {
      return 1;
    }
    if (ln < 0)
    {

      return -1;
    }
    return 0;
  }
  int Orientation3D(Cartesian p1, Cartesian p2)
  {
    Cartesian a = Sub(p1, this);
    Cartesian b = Sub(p2, this);
    Cartesian sa = a.CrossProduct(a, b);
    double ln = sa.x + sa.y + sa.z;
    if (ln > 0)
    {
      return 1;
    }
    if (ln < 0)
    {

      return -1;
    }
    return 0;
  }
  boolean Parallel2D(Cartesian p)
  {
    Cartesian na = Rotate90Degrees2D();
    return Math.abs(DotProduct2D(na, p)) < Epsilon;
  }
  boolean Parallel3D(Cartesian p)
  {
    Cartesian na = CrossProduct(this, p);
    return Math.abs(DotProduct(na, p)) < Epsilon;
  }
  boolean Parallel2D(Cartesian p, double r)
  {
    Cartesian na = Rotate90Degrees2D();
    return r * r < Math.abs(DotProduct2D(na, p));
  }
  Cartesian Rotate90Degrees2D()
  {
    return new Cartesian(-y, x);
  }
  Cartesian Flip2D()
  {
    return new Cartesian(x, -y);
  }
  Cartesian Mirror2D()
  {
    return new Cartesian(-x, y);
  }

  boolean Collides2D(Cartesian p)
  {
    return Math.abs(x - p.x) < Epsilon && Math.abs(y - p.y) < Epsilon;
  }
  boolean Collides2D(Cartesian p, double threshold)
  {
    return Math.abs(x - p.x) < threshold && Math.abs(y - p.y) < threshold;
  }
  boolean Collides3D(Cartesian p)
  {
    return Math.abs(x - p.x) < Epsilon && Math.abs(y - p.y) < Epsilon && Math.abs(z - p.z) < Epsilon;
  }
  boolean Collides3D(Cartesian p, double threshold)
  {
    return Math.abs(x - p.x) < threshold && Math.abs(y - p.y) < threshold && Math.abs(z - p.z) < threshold;
  }
  Normal Normal()
  {
    double len = Length();
    double im = 1.0f;
    if (len != 0) im /= len; 
    else return new Normal(0, 0, 0);
    return new Normal(x * im, y * im, z * im);
  }
  Vector UnityVector2D()
  {
    double len = Length2D();
    double im = 1.0f;
    if (len != 0) im /= len;
    return new Vector(x * im, y * im);
  }
  Normal Normal2D()
  {
    double len = Length2D();
    double im = 1.0f;
    if (len != 0) im /= len;
    return new Normal(x * im, y * im);
  }
  Cartesian Rotate2D(double ang)
  {
    double cs = (double)Math.cos(ang);
    double sn = (double)Math.sin(ang);
    return new Cartesian(x * cs - y * sn, y * cs + x * sn);
  }
  Cartesian Rotate3D(double phi, double theta)
  {
    double r = (double)Math.sqrt(x * x + y * y);
    double p = Length();
    theta += (double)Math.atan2(y, x);
    phi += (double)Math.acos(z / p);
    return new Cartesian(p * (double)Math.sin(phi) * (double)Math.cos(theta), p * (double)Math.sin(phi) * (double)Math.sin(theta), p * (double)Math.cos(phi));
  }
  Cartesian Rotate2D(Cartesian at, double ang)
  {
    double cs = (double)Math.cos(ang);
    double sn = (double)Math.sin(ang);
    return Add(at, new Cartesian(x * cs - y * sn, y * cs + x * sn) );
  }
  Cartesian Rotate3D(Cartesian at, double phi, double theta)
  {
    double r = (double)Math.sqrt(x * x + y * y);
    double p = Length();
    theta += (double)Math.atan2(y, x);
    phi += (double)Math.acos(z / p);
    return Add(at, new Cartesian(p * (double)Math.sin(phi) * (double)Math.cos(theta), p * (double)Math.sin(phi) * (double)Math.sin(theta), p * (double)Math.cos(phi)) );
  }
  Cartesian Unity()
  {
    double len = Length();
    double im = 1.0f;
    if (len != 0) im /= len;
    return new Cartesian(x * im, y * im, z * im);
  }
  Cartesian Unity2D()
  {
    double len = Length2D();
    double im = 1.0f;
    if (len != 0) im /= len;
    return new Cartesian(x * im, y * im);
  }
  Normal Normal3D()
  {
    double len = Length();
    double im = 1.0f;
    if (len != 0) im /= len;
    return new Normal(x * im, y * im, z * im);
  }
  void Normalize()
  {
    double len = Length();
    double im = 1.0f;
    if (len != 0) im /= len;
    x = (x * im);
    y = y * im;
    z = z * im;
  }
  boolean EqualCube(Cartesian p, double threshold)
  {
    return Math.abs(x - p.x) < threshold && Math.abs(y - p.y) < threshold && Math.abs(z - p.z) < threshold;
  }
  boolean Equals(Cartesian p0, Cartesian p1)
  {
    return p0 == p1;
  }
  double EnclosedAngle2D(Cartesian p0, Cartesian p1)
  {
    Cartesian u1 = p0.Unity2D();
    Cartesian u2 = p1.Unity2D();
    double dp = DotProduct2D(u1, u2);
    return (double)Math.acos(dp);
  }
  double EnclosedAngle(Cartesian p0, Cartesian p1)
  {
    Cartesian u1 = p0.Unity();
    Cartesian u2 = p1.Unity();
    double dp = DotProduct(u1, u2);
    return (double)Math.acos(dp);
  }
  boolean Collision(double r1, Cartesian p, double r2)
  {
    Cartesian a = Sub(this, p);
    double r3 = r1 * r1 + r2 * r2;
    if (a.LengthSquared() < r3) return true;
    return false;
  }
  boolean Collision(Cartesian p1, double r1, Cartesian p2, double r2)
  {
    Cartesian a = Sub(p1, p2);
    double r3 = r1 * r1 + r2 * r2;
    if (a.LengthSquared() < r3) return true;
    return false;
  }
  byte ClassifyPosition2D(Cartesian p0, Cartesian p1)
  {
    Cartesian p2 = this;
    Cartesian a = Sub(p1, p0);
    Cartesian b = Sub(p2, p0);
    double sa = a.x * b.y - a.y * b.x;
    if (sa > Epsilon) return PositionLeft;
    if (sa < -Epsilon) return PositionRight;
    if ((a.x * b.x) < -Epsilon || (a.y * b.y) < -Epsilon) return PositionBehind;
    if (a.LengthSquared() < b.LengthSquared())
    {
      return PositionBeyond;
    }
    if (p0.Equal(p0, p2))
    {
      return PositionOrigin;
    }
    if (p1.Equal(p1, p2))
    {
      return PositionDestination;
    }
    return PositionBetween;
  }
  byte ClassifyPosition2D(Point p0, Point p1)
  {
    Point p2 = AsPoint(1);
    Vector a = p1.Sub(p1, p0);
    Vector b = p2.Sub(p2, p0);
    double sa = a.x * b.y - a.y * b.x;
    if (sa > Epsilon) return PositionLeft;
    if (sa < -Epsilon) return PositionRight;
    if ((a.x * b.x) < -Epsilon || (a.y * b.y) < -Epsilon) return PositionBehind;
    if (a.LengthSquared() < b.LengthSquared())
    {
      return PositionBeyond;
    }
    if (Equal(p0, p2))
    {
      return PositionOrigin;
    }
    if (Equal(p1, p2))
    {
      return PositionDestination;
    }
    return PositionBetween;
  }

  // Geometry Internal Functions
  Cartesian Add(Cartesian v, Vector p)
  {
    return new Cartesian(v.x + p.x, v.y + p.y, v.z + p.z);
  }
  Cartesian Sub(Cartesian v, Vector p)
  {
    return new Cartesian(v.x - p.x, v.y - p.y, v.z - p.z);
  }
  Cartesian Add(Cartesian v, Normal p)
  {
    return new Cartesian(v.x + p.x, v.y + p.y, v.z + p.z);
  }
  Cartesian Sub(Cartesian v, Normal p)
  {
    return new Cartesian(v.x - p.x, v.y - p.y, v.z - p.z);
  }
  Cartesian Add(Cartesian v, Cartesian p)
  {
    return new Cartesian(v.x + p.x, v.y + p.y, v.z + p.z);
  }
  Cartesian Sub(Cartesian v, Cartesian p)
  {
    return new Cartesian(v.x - p.x, v.y - p.y, v.z - p.z);
  }
  Cartesian Sub(Cartesian a, Point p)
  {
    return new Cartesian(a.x - p.x, a.y - p.y, a.z - p.z);
  }
  Cartesian Negative(Cartesian p)
  {
    return new Cartesian(-p.x, -p.y, -p.z);
  }
  double Distance(Cartesian p1, Cartesian p2)
  {
    return (Sub(p1, p2)).Length();
  }
  double DistanceSquared(Cartesian p1, Cartesian p2)
  {
    return (Sub(p1, p2)).LengthSquared();
  }
  boolean Equal(Cartesian a, Cartesian p)
  {
    if (a==(null) && p==(null)) { 
      return true;
    } else if (a==(null) || p==(null)) { 
      return false;
    } else
      return Math.abs(a.x - p.x) < Epsilon && Math.abs(a.y - p.y) < Epsilon && Math.abs(a.z - p.z) < Epsilon;
  }
  boolean NotEqual(Cartesian a, Cartesian p)
  {
    return !Equal(a, p);
  }
  boolean LessThan(Cartesian a, Cartesian p)
  {
    return (a.x < p.x) || (a.x == p.x && a.y < p.y);
  }
  boolean Greater(Cartesian a, Cartesian p)
  {
    return (a.x > p.x) || (a.x == p.x && a.y > p.y);
  }
  boolean Parallel3D(Cartesian p, double r)
  {
    Cartesian na = CrossProduct(this, p);
    return r * r < Math.abs(DotProduct(p, na));
  }
  double Sine(double invr)
  {
    return x * invr;
  }
  double Cosine(double invr)
  {
    return y * invr;
  }
  double LengthSquared()
  {
    return x * x + y * y + z * z;
  }
  double Length()
  {
    return (double)Math.sqrt(LengthSquared());
  }
  double InvLength()
  {
    return (double)(1.0f / Math.sqrt(LengthSquared()));
  }
  double Tangent(byte axis)
  {
    double t = 0;
    if (axis == AxisZ)
    {
      if (x != 0)
      {
        t = (double)y / x;
      }
    } else if (axis == AxisY)
    {
      if (x != 0)
      {
        t = (double)z / x;
      }
    } else
    {
      if (z != 0)
      {
        t = (double)y / z;
      }
    }
    return t;
  }

  double PolarAngle2D()
  {
    if ((x == 0.0) && (y == 0.0))
    {
      return -1;
    }
    if (x == 0.0)
    {
      return ((y > 0) ? PI90 : PI270);
    }
    double theta = Math.atan(y / x);
    if (x > 0)
    {
      return ((y >= 0) ? theta : PI2 + theta);
    } else
    {
      return (Math.PI + theta);
    }
  }
  Point AsPoint(double r)
  {
    return new Point(x, y, z, r);
  }
  Normal AsNormal()
  {
    return new Normal(x, y, z);
  }
  Vector AsVector()
  {
    return new Vector(x, y, z);
  }
  double Angle()
  {
    double theta = (double)(Math.atan2(y, x));
    return theta;
  }
  double Length2DSquared()
  {
    return x * x + y * y;
  }
  double Length2D()
  {
    return (double)Math.sqrt(Length2DSquared());
  }
  double Dot(Cartesian n1, Cartesian n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double DotProduct(Cartesian n1, Cartesian n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double AbsDot(Cartesian n1, Cartesian n2)
  {
    return Math.abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
  }
  double Dot(Cartesian n1, Point n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double DotProduct(Cartesian n1, Point n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double AbsDot(Cartesian n1, Point n2)
  {
    return Math.abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
  }
  double Dot(Cartesian n1, Vector n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double DotProduct(Cartesian n1, Vector n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double AbsDot(Cartesian n1, Vector n2)
  {
    return Math.abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
  }
  double Dot(Cartesian n1, Normal n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double DotProduct(Cartesian n1, Normal n2)
  {
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  }
  double AbsDot(Cartesian n1, Normal n2)
  {
    return Math.abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
  }
  void Scale2D(double s)
  {
    x *= s;
    y *= s;
  }
  void Scale(double s)
  {
    x *= s;
    y *= s;
    z *= s;
  }
  Vector Scale(Cartesian a, double s) {
    return new Vector(a.x * s, a.y * s, a.z * s);
  }
  Vector Project2D(Cartesian onto)
  {
    double dp = DotProduct2D(onto, onto);
    if (dp > 0)
    {
      double pd = DotProduct2D(this, onto);
      return onto.Scaled(pd / dp).AsVector();
    }
    return new Vector(onto.x, onto.y, onto.z);
  }
  Vector Project(Cartesian onto)
  {
    double dp = DotProduct(onto, onto);
    if (dp > 0)
    {

      double pd = DotProduct(this, onto);
      return onto.Scaled(pd / dp).AsVector();
    }
    return new Vector(onto.x, onto.y, onto.z);
  }
  //
  Vector Project2D(Cartesian p0, Cartesian onto)
  {
    double dp = DotProduct2D(onto, onto);
    if (dp > 0)
    {
      double pd = DotProduct2D(p0, onto);
      return onto.Scaled(pd / dp).AsVector();
    }
    return new Vector(onto.x, onto.y, onto.z);
  }
  Vector Project(Cartesian p0, Cartesian onto)
  {
    double dp = DotProduct(onto, onto);
    if (dp > 0)
    {

      double pd = DotProduct(p0, onto);
      return onto.Scaled(pd / dp).AsVector();
    }
    return new Vector(onto.x, onto.y, onto.z);
  }
  Cartesian Scaled(Cartesian v) {

    return new Vector( x*v.x, y*v.y, z*v.z);
  }
  Cartesian Scaled(double v) {
    return new Cartesian(x*v, y*v, z*v);
  }
  double DotProduct2D(Cartesian p0, Cartesian p1)
  {
    return (p0.x * p1.x + p0.y * p1.y);
  }
  double DotProduct2DXY(Cartesian p0, Cartesian p1)
  {
    return (p0.x * p1.x + p0.y * p1.y);
  }
  double DotProduct2DYZ(Cartesian p0, Cartesian p1)
  {
    return (p0.z * p1.z + p0.y * p1.y);
  }
  double DotProduct2DZX(Cartesian p0, Cartesian p1)
  {
    return (p0.x * p1.x + p0.z * p1.z);
  }
  double ScalarProduct(Cartesian p0, Cartesian p1)
  {
    return (p0.x * p1.x + p0.y * p1.y + p0.z * p1.z);
  }
  Cartesian Cross(Cartesian v1, Cartesian v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian CrossProduct(Cartesian v1, Cartesian v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian VectorProduct(Cartesian v1, Cartesian v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }

  Cartesian Cross(Cartesian v1, Normal v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian VectorProduct(Cartesian v1, Normal v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian CrossProduct(Cartesian v1, Normal v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian Cross(Cartesian v1, Vector v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian VectorProduct(Cartesian v1, Vector v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian CrossProduct(Cartesian v1, Vector v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian Cross(Cartesian v1, Point v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian VectorProduct(Cartesian v1, Point v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian CrossProduct(Cartesian v1, Point v2)
  {
    return new Cartesian((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
  }
  Cartesian Negative() {
    return new Cartesian(-x, -y, -z);
  }
}
class Vector extends Cartesian {
  double w;
  Vector() {
    w=1;
  }
  Vector(double xx, double yy) {
    x=xx;
    y=yy;
    z=0;
    w=1;
  }
  Vector(double xx, double yy, double zz) {
    x=xx;
    y=yy;
    z=zz;
    w=1;
  }
  Vector(double xx, double yy, double zz, double ww) {
    x=xx;
    y=yy;
    z=zz;
    w=ww;
  }
  Vector(Vector v) {
    x=v.x;
    y=v.y;
    z=v.z;
    w=v.w;
  }
  double Length() {
    return Math.sqrt(x*x+y*y+z*z);
  }
  double InvLength() {
    return 1.0/Math.sqrt(x*x+y*y+z*z);
  }
  double LengthSqrd() {
    return (x*x+y*y+z*z);
  }
  Normal Normal() {
    double vi=InvLength();
    return new Normal(vi*x, vi*y, vi*z);
  }
  Vector Negative() {
    return new Vector(-x, -y, -z);
  }
  Vector Absolute() {
    return new Vector(Math.abs(x), Math.abs(y), Math.abs(z));
  }
  void Normalize() {
    double vi=InvLength();
    Scale(vi);
  }
  void Add(Vector v) {
    x+=v.x;
    y+=v.y;
    z+=v.z;
  }
  void Sub(Vector v) {
    x-=v.x;
    y-=v.y;
    z-=v.z;
  }
  void Scale(Vector v) {
    x*=v.x;
    y*=v.y;
    z*=v.z;
  }
  void Scale(double v) {
    x*=v;
    y*=v;
    z*=v;
  }

  Vector Scaled(Vector v) {

    return new Vector( x*v.x, y*v.y, z*v.z);
  }
  Vector Scaled(double v) {
    return new Vector( x*v, y*v, z*v);
  }
  void Clamp(double from, double to) {
    if (x<from) {
      x=from;
    } else if (x>to) {
      x=to;
    }
    if (y<from) {
      y=from;
    } else if (y>to) {
      y=to;
    }
    if (z<from) {
      z=from;
    } else if (z>to) {
      z=to;
    }
  }
  double Dot(Vector v) {
    return x*v.x+y*v.y+z*v.z;
  }
  Vector Cross(Vector v) {
    return new Vector(y*v.z+z*v.y, x*v.z+z*v.x, y*v.x+y*v.x);
  }
}
class Matrix {
  double m[][]=new double[4][4];
  Matrix() {
    m[0][0]=1;
    m[1][1]=1;
    m[2][2]=1;
    m[3][3]=1;
  }
  Matrix RotateX(double a) {
    Matrix x=new Matrix();
    x.m[1][1]=Math.cos(a);
    x.m[1][2]=Math.sin(a);
    x.m[2][1]=-Math.sin(a);
    x.m[2][2]=Math.cos(a);
    return x;
  }
  Matrix RotateY(double a) {
    Matrix y=new Matrix();
    y.m[0][0]=Math.cos(a);
    y.m[2][0]=-Math.sin(a);
    y.m[0][2]=Math.sin(a);
    y.m[2][2]=Math.cos(a);
    return y;
  }
  Matrix RotateZ(double a) {

    Matrix z=new Matrix();
    z.m[0][0]=Math.cos(a);
    z.m[1][0]=Math.sin(a);
    z.m[0][1]=-Math.sin(a);
    z.m[1][1]=Math.cos(a);
    return z;
  }
  Matrix Mult(Matrix mult, Matrix by) {

    Matrix matrix=new Matrix();
    for (int j=0; j<4; j++) {
      for (int i=0; i<4; i++) {
        matrix.m[i][j]=mult.m[i][0]*by.m[0][j]+mult.m[i][1]*by.m[1][j]+mult.m[i][2]*by.m[2][j]+mult.m[i][3]*by.m[3][j];
      }
    }
    return matrix;
  }
  Vector Mult(Matrix mult, Vector by) {
    double x=by.x*mult.m[0][0]+by.y*mult.m[1][0]+by.z*mult.m[2][0]+by.w*mult.m[3][0];
    double y=by.x*mult.m[0][1]+by.y*mult.m[1][1]+by.z*mult.m[2][1]+by.w*mult.m[3][1];
    double z=by.x*mult.m[0][2]+by.y*mult.m[1][2]+by.z*mult.m[2][2]+by.w*mult.m[3][2];
    double w=by.x*mult.m[0][3]+by.y*mult.m[1][3]+by.z*mult.m[2][3]+by.w*mult.m[3][3];
    return new Vector(x, y, z, w);
  }
}

    class Quarternion {
        double[] q = new double[4];
        Quarternion()
        {
            q[0] = 1; q[1] = 1; q[2] = 1;
            q[3] = 0;
        }
        Quarternion(Quarternion copy)
        {
            q[0] = copy.q[0];
            q[1] = copy.q[1];
            q[2] = copy.q[2];
            q[3] = copy.q[3];
        }
         Quarternion(double qw, double qx, double qy, double qz)
        {
            q[0] = qx; q[1] = qy; q[2] = qz; q[3] = qw;
        }
        Quarternion(Vector v, double qs)
        {
            q[0] = v.x; q[1] = v.y; q[2] = v.z; q[3] = qs;
        }
        Quarternion(Vector v)
        {
            q[0] = v.x; q[1] = v.y; q[2] = v.z; q[3] = 0;
        }
        Quarternion(Normal v)
        {
            q[0] = v.x; q[1] = v.y; q[2] = v.z; q[3] = 0;
        }

        Quarternion(Point v)
        {
            q[0] = v.x; q[1] = v.y; q[2] = v.z; q[3] = 0;
        }
        Quarternion(double a, Normal v)
        {
            double a2 = a * 0.5;
            double sa2 = Math.sin(a2);
            double ca2 = Math.cos(a2);
            q[0] = v.x * sa2; q[1] = v.y * sa2; q[2] = v.z * sa2; q[3] = ca2;
        }
        Quarternion(double a, Point v)
        {
            double a2 = a * 0.5;
            double sa2 = Math.sin(a2);
            double ca2 = Math.cos(a2);
            q[0] = v.x * sa2; q[1] = v.y * sa2; q[2] = v.z * sa2; q[3] = ca2;
        }
        Quarternion Mult(Quarternion a, Quarternion b) {
            Vector q = new Vector(a.q[0], a.q[1], a.q[2]);
            Vector p = new Vector(b.q[0], b.q[1], b.q[2]);
            Vector pqx = q.Scaled( b.q[3]);
            Vector pqy = p.Scaled(a.q[3]); 
            Vector pqz = p.Cross(p, q).AsVector();
            Vector pq=pqx.Add(pqx,pqx.Add(pqy,pqz)).AsVector();
            double pqs = b.q[3] * a.q[3] - p.DotProduct(p, q);
            return new Quarternion(pqs,pq.x,pq.y,pq.z);
        }
        Quarternion Add(Quarternion a, Quarternion b)
        {
            return new Quarternion(a.q[0] + b.q[0], a.q[1] + b.q[1], a.q[2] + b.q[2], a.q[3] + b.q[3]);
        }
        Quarternion Sub(Quarternion a, Quarternion b)
        {
            return new Quarternion(a.q[0] - b.q[0], a.q[1] - b.q[1], a.q[2] - b.q[2], a.q[3] - b.q[3]);
        }
        Quarternion Scale(Quarternion a, double s)
        {
            Quarternion sc = new Quarternion(a);
            sc.q[0] *= s;
            sc.q[1] *= s;
            sc.q[2] *= s;
            sc.q[3] *= s;
            return sc;
        }
        Quarternion Scaled(double s)
        {
            Quarternion sc = new Quarternion(this);
            sc.q[0] *= s;
            sc.q[1] *= s;
            sc.q[2] *= s;
            sc.q[3] *= s;
            return sc;
        }
        Quarternion Perspect(Quarternion a, double d)
        {
            double s = 1.0 / d;
            Quarternion sc = new Quarternion(a);
            sc.q[0] *= s;
            sc.q[1] *= s;
            sc.q[2] *= s;
            sc.q[3] *= s;
            return sc;
        }
        Vector Rotate(Vector v)
        {
            Quarternion vq = new Quarternion(v);
            Quarternion p = this.Mult((this.Mult(this, vq)), InverseOfUnit());
            return new Vector(p.q[0], p.q[1], p.q[2], p.q[3]);
        }
        Normal Rotate(Normal v)
        {
            Quarternion vq = new Quarternion(v);
            Quarternion p = this.Mult((this.Mult(this, vq)), InverseOfUnit());
            return new Normal(p.q[0], p.q[1], p.q[2]);
        }
        Point Rotate(Point v)
        {
            Quarternion vq = new Quarternion(v);
            Quarternion p = this.Mult((this.Mult(this, vq)), InverseOfUnit());
            return new Point(p.q[0], p.q[1], p.q[2]);
        }
        void Normalize()
        {
            Quarternion c = Conjugate();
            double sq = 1.0 / Math.sqrt(q[0] * c.q[0] + q[1] * c.q[1] + q[2] * c.q[2] + q[3] * c.q[3]);
            q[0] *= sq;
            q[1] *= sq;
            q[2] *= sq;
            q[3] *= sq;
        }
        double MagnitudeSquared()
        {
            double sq = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            return sq;
        }
        double Magnitude()
        {
            double sq = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            return sq;
        }
        double InvMagnitude()
        {
            double sq = 1.0 / Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            return sq;
        }
        double ScaleMagnitude()
        {
            double sq = 2.0 / Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            return sq;
        }
        double ScaleMagnitude(double s)
        {
            double sq = s / Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
            return sq;
        }
        Quarternion Normal() {
            Quarternion n = new Quarternion(this);
            n.Normalize();
            return n;
        }
        Quarternion Conjugate() {
            return new Quarternion(q[3], -q[0], -q[1], -q[2]);
        }
        Quarternion InverseOfUnit()
        {
            return new Quarternion(q[3], -q[0], -q[1], -q[2]);
        }
        Quarternion Inverse()
        {
            Quarternion qn = new Quarternion(this);
            qn.Normalize();
            double ns = 1.0 / this.MagnitudeSquared();
            return qn.Conjugate().Scaled(ns);
        }
        Quarternion Net(Quarternion sq, Quarternion[] qs)
        {
            Quarternion e = new Quarternion(sq);
            for(int i=0;i<qs.length;i++){
            Quarternion f =qs[i];
              e=e.Mult(e,f);
            }
            return e;
        }
        Quarternion LinearInterop(Quarternion a, Quarternion b, double t)
        {  Quarternion ar=a.Scaled(1.0 - t);
           Quarternion br=b.Scaled(t);
            Quarternion r = (ar.Add(ar,br));
            r.Normalize();
            return r;
        }
        Quarternion SphereInterop(Quarternion a, Quarternion b, double t)
        {
            Point q = new Point(a.q[0], a.q[1], a.q[2]);
            Point p = new Point(b.q[0], b.q[1], b.q[2]);
            //
            double phi = Math.acos(p.DotProduct(q, p) * q.InvLength());
            double sp= Math.sin(phi);
            double wa = (Math.sin(1.0 - t) * phi) / sp;
            double wb = Math.sin(t * phi) /sp ;
            Quarternion ar=a.Scaled(wa);
            Quarternion br=b.Scaled(wb);
            Quarternion r = ar.Add(ar, br);
            //r.Normalize();
            return r;
        }
        Matrix GetMatrix()
        {
            Matrix qm = new Matrix();
            qm.m[0][ 0] = 1.0 - 2.0 * (q[1] * q[1] - q[2] * q[2]); qm.m[1][ 0] = 2.0 * (q[0] * q[1] + q[2] * q[3]); qm.m[2][ 0] = 2.0 * (q[0] * q[2] - q[1] * q[3]);
            //
            qm.m[0][ 1] = 2.0 * (q[0] * q[1] - q[2] * q[3]); qm.m[1][ 1] = 1.0 - 2.0 * (q[0] * q[0] + q[2] * q[2]); qm.m[2][ 1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
            //
            qm.m[0][ 2] = 2.0 * (q[0] * q[2] + q[1] * q[3]); qm.m[1][ 2] = 2.0 * (q[1] * q[2] - q[0] * q[3]); qm.m[2][2] = 1.0 - 2.0 * (q[0] * q[0] - q[1] * q[1]);
            return qm;
        }
         Matrix GetMatrix(double scale)
        {
            Matrix qm = new Matrix();
            qm.m[0][ 0] = 1.0 - scale * (q[1] * q[1] - q[2] * q[2]); qm.m[1][ 0] = scale * (q[0] * q[1] + q[2] * q[3]); qm.m[2][ 0] = scale * (q[0] * q[2] - q[1] * q[3]);
            //
            qm.m[0][ 1] = scale * (q[0] * q[1] - q[2] * q[3]); qm.m[1][ 1] = 1.0 - scale * (q[0] * q[0] + q[2] * q[2]); qm.m[2][ 1] = scale * (q[1] * q[2] + q[0] * q[3]);
            //
            qm.m[0][ 2] = scale * (q[0] * q[2] + q[1] * q[3]); qm.m[1][ 2] = scale * (q[1] * q[2] - q[0] * q[3]); qm.m[2][ 2] = 1.0 - scale * (q[0] * q[0] - q[1] * q[1]);

            return qm;
        }
        Quarternion(Matrix m)
        {
            double trace = m.m[0][ 0] + m.m[1][ 1] + m.m[2][ 2];
            if (trace > 0) {
                double s = Math.sqrt(trace + 1.0);
                q[3] = s * 0.5;
                double t = 0.5 / s;
                q[0] = (m.m[2][ 1] - m.m[1][ 2]) * t;
                q[1] = (m.m[0][ 2] - m.m[2][ 0]) * t;
                q[2] = (m.m[1][ 0] - m.m[0][ 1]) * t;
            } else {
                int i = 0;
                if (m.m[1][ 1] > m.m[0][ 0]) { i = 1; }
                if (m.m[2][ 2] > m.m[i][ i]) { i = 2; }
                int[] next = new int[] { 1, 2, 0 };
                int j = next[i];
                int k = next[j];

                double s = Math.sqrt((m.m[i][ j] - (m.m[j][ j] + m.m[k][ k]))+1.0);
                q[i] = s * 0.5;
                double t;
                if (s != 0) { t = 0.5 / s; } else { t = s; }
                q[3] = (m.m[k][ j] - m.m[j][ k]) * t;
                q[j] = (m.m[j][ i] + m.m[i][ j]) * t;
                q[k] = (m.m[k][ i] + m.m[i][ k]) * t;

            }
        }
    }
    class Anim//sqt
    {
        Vector scale;
        Quarternion transform;
        Vector translate;
    }
    class AxisAngle//sqt
    {
        Vector angles;
        Vector direction;
        Quarternion transform;
    }
    class Orientation
    {
        Vector angle;
        Normal direction;
        Point at;
        Vector scale;
        Vector translate;
        Quarternion transform;
    }
    public class Dual
    {
        Quarternion a;
        Quarternion b;
        Quarternion GetValue()
        {
            return a.Add(a,b);
        }
    }
    public class SphereTransform
    {
        Quarternion x;
        Quarternion y;
        Quarternion z;
        Quarternion GetTransform(Vector a)
        {
            x = new Quarternion(a.x, new Normal(1, 0, 0));
            y = new Quarternion(a.y, new Normal(0, 1, 0));
            z = new Quarternion(a.z, new Normal(0, 0, 1));
            return x.Mult(x,x.Mult(y, z));
        }
        Quarternion GetTransform(Vector a, SphereTransform to, double t)
        {

            x = new Quarternion(a.x, new Normal(1, 0, 0));
            y = new Quarternion(a.y, new Normal(0, 1, 0));
            z = new Quarternion(a.z, new Normal(0, 0, 1));
            return x.SphereInterop(x.Mult(x,x.Mult(y, z)), to.GetTransform(),t);
        }
        Quarternion GetTransform(Vector a, Vector b, SphereTransform to, double t)
        {

            x = new Quarternion(a.x, new Normal(1, 0, 0));
            y = new Quarternion(a.y, new Normal(0, 1, 0));
            z = new Quarternion(a.z, new Normal(0, 0, 1));
            return x.SphereInterop(x.Mult(x,x.Mult(y, z)), to.GetTransform(b), t);
        }
        Quarternion GetTransform(Vector a, Quarternion to, double t)
        {

            x = new Quarternion(a.x, new Normal(1, 0, 0));
            y = new Quarternion(a.y, new Normal(0, 1, 0));
            z = new Quarternion(a.z, new Normal(0, 0, 1));
            return x.SphereInterop(x.Mult(x,x.Mult(y, z)), to, t);
        }
        Quarternion GetTransform()
        {
            return x.Mult(x,x.Mult(y, z));
        }

    }
