class Ray
{
  public:
    Ray(const Vector &NewC, const Vector &Newu)
    {
      C = NewC;
      u = Newu;
    }

    Vector C, u;
};
