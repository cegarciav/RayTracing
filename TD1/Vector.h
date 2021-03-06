class Vector
{
  public:
    Vector(double x=0, double y=0, double z=0)
    {
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
    }

    Vector operator+(const Vector &B)
    {
      return Vector(coords[0] + B[0],
                    coords[1] + B[1],
                    coords[2] + B[2]
                    );
    }

    Vector operator-(const Vector &B)
    {
      return Vector(coords[0] - B[0],
                    coords[1] - B[1],
                    coords[2] - B[2]
                    );
    }

    friend Vector operator*(const double a, const Vector &B)
    {
      return Vector(a * B[0], a * B[1], a * B[2]);
    }

    Vector &operator=(const Vector &B)
    {
        coords[0] = B[0];
        coords[1] = B[1];
        coords[2] = B[2];
        return *this;
    }

    double operator[](int cord) const
    {
      return coords[cord];
    }

    double operator[](int cord)
    {
      return coords[cord];
    }

    double getNorm2()
    {
      return coords[0] * coords[0] +
              coords[1] * coords[1] +
              coords[2] * coords[2];
    }

    double dot(const Vector& B)
    {
      return coords[0] * B[0] + coords[1] * B[1] + coords[2] * B[2];
    }

    Vector getNormalized()
    {
      Vector n_vector = *this;
      n_vector = (1 / sqrt(n_vector.getNorm2())) * n_vector;
      return n_vector;
    }

    double coords[3];

};
