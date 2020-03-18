class Ray
{
  public:
  	Ray(){};
    Ray(const Vector &NewC, const Vector &Newu)
    {
      C = NewC;
      u = Newu;
    }

    Vector C, u;
};
