class Triangle : public Object
{
  public:
    Triangle(
              const Vector &PointA,
              const Vector &PointB,
              const Vector &PointC,
              const Vector &NewAlbedo=Vector(0.,0.,0.),
              bool mirror=false,
              bool transparent=false,
              int normal_sign = 1
            )
    {
      A = PointA;
      B = PointB;
      C = PointC;
      N = normal_sign * (B - A).cross(C - A).getNormalized();
      albedo = NewAlbedo;
      is_mirror = mirror;
      is_transparent = transparent;
    }

    bool intersect(const Ray &r, Vector &Point, Vector &Normal, double &t) const
    {
      double alpha, beta, gamma;
      return intersect(r, Point, Normal, t, alpha, beta, gamma);
    }

    bool intersect(
                    const Ray &r, Vector &Point, Vector &Normal,
                    double &t, double &alpha, double &beta,
                    double &gamma
                  ) const
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

      // Compute vectors        
      Vector AC = C - A,
              AB = B - A,
              AP = Point - A;

      // Compute dot products
      double AC2 = AC.getNorm2(),
              AC_AB = AC.dot(AB),
              AC_AP = AC.dot(AP),
              AB2 = AB.getNorm2(),
              AB_AP = AB.dot(AP);

      // Compute barycentric coordinates
      double inv_det_M = 1 / (AB2 * AC2 - AC_AB * AC_AB);
      beta = (AB_AP * AC2 - AC_AB * AC_AP) * inv_det_M;
      gamma = (AB2 * AC_AP - AB_AP * AC_AB) * inv_det_M;
      alpha = 1 - beta - gamma;

      // Check if point is in triangle
      return (beta >= 0) && (gamma >= 0) && (beta + gamma < 1);
    }

    Vector A, B, C, N;
};