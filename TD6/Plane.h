class Plane : public Object
{
  public:
    Plane(
              const Vector &PointA,
              const Vector &Normal,
              const Vector &NewAlbedo=Vector(0.,0.,0.),
              bool mirror=false,
              bool transparent=false,
              int normal_sign = 1
            )
    {
      A = PointA;
      N = Normal;
      albedo = NewAlbedo;
      is_mirror = mirror;
      is_transparent = transparent;
    }

    bool intersect(const Ray &r, Vector &Point, Vector &Normal, double &t, Vector &colour) const
    {
      double b = r.u.dot(N);
      if (b == 0 || isnan(b))
        return false;
      double a = (C - r.C).dot(N);
      t = a / b;
      if (t < 0)
        return false;

      Point = r.C + t * r.u;
      Normal = N;
      colour = albedo;

      return true;
    }

    Vector A, B, C, N;
};