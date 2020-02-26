class Sphere
{
  public:
    Sphere(const Vector &NewO, double NewR, const Vector &NewAlbedo)
    {
      O = NewO;
      R = NewR;
      albedo = NewAlbedo;
    }

    bool intersect(Ray &r, Vector &Point, Vector &Normal, double &t)
    {
      // solves at^2 + bt + c
      double a = 1;
      double b = 2 * r.u.dot(r.C - O);
      double c  = (r.C - O).getNorm2() - R * R;

      //determinant < 0
      double delta = b * b - 4 * a * c;
      if (delta < 0)
        return false;

      //else, solve the equation
      double t1 = (-b + sqrt(delta)) / (2 * a);
      double t0 = (-b - sqrt(delta)) / (2 * a);

      //if t1 (the furthest) is negative, t0 is negative too
      if (t1 < 0)
        return false;

      if (t0 > 0)
        t = t0;
      else
        t = t1;

      Point = r.C + t * r.u;
      Normal = (Point - O).getNormalized();

      return true;
    }

    double getAlbedo(int index)
    {
      return albedo[index];
    }

    Vector O;
    double R;
    Vector albedo;
};
