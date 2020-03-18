class BBox
{
  public:
    BBox(){};
    BBox(const Vector &nbmin, const Vector &nbmax)
    {
      bmin = nbmax;
      bmax = nbmax;
    }

    bool intersect(const Ray &r) const
    {
      double t_1x = (bmin[0] - r.C[0]) / r.u[0];
      double t_2x = (bmax[0] - r.C[0]) / r.u[0];
      double t_xmin = std::min(t_1x, t_2x);
      double t_xmax = std::max(t_1x, t_2x);

      double t_1y = (bmin[1] - r.C[1]) / r.u[1];
      double t_2y = (bmax[1] - r.C[1]) / r.u[1];
      double t_ymin = std::min(t_1y, t_2y);
      double t_ymax = std::max(t_1y, t_2y);

      double t_1z = (bmin[2] - r.C[2]) / r.u[2];
      double t_2z = (bmax[2] - r.C[2]) / r.u[2];
      double t_zmin = std::min(t_1z, t_2z);
      double t_zmax = std::max(t_1z, t_2z);

      return (
                std::min(std::min(t_xmax, t_ymax), t_zmax) > 0 &&
                std::min(std::min(t_xmax, t_ymax), t_zmax) - 
                std::max(std::max(t_xmin, t_ymin), t_zmin) > 0
              );
    }

    Vector bmin, bmax;
};
