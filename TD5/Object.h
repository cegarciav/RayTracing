class Object
{
  public:
    Object(){};
    virtual bool intersect(const Ray &r, Vector &Point, Vector &Normal, double &t) const = 0;

    Vector albedo;
    bool is_mirror;
    bool is_transparent;
    bool is_light;
};
