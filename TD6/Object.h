class Object
{
  public:
    Object(){};
    virtual bool intersect(const Ray &r, Vector &Point, Vector &Normal, double &t, Vector &colour) const = 0;

    Vector albedo;
    double phong_exp;
    double ks;
    bool is_mirror;
    bool is_transparent;
    bool is_light;
};
