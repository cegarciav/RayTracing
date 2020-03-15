class Scene
{
  public:
    Scene(const Sphere *light, double new_intensity)
    {
      lumiere = light;
      intensity = new_intensity;
    }
    
    void addSphere(const Sphere &s)
    {
      objects.push_back(&s);
    }

    void addTriangle(const Triangle &t)
    {
      objects.push_back(&t);
    }

    void addGeometry(const Geometry &g)
    {
      objects.push_back(&g);
    }

    bool intersect(const Ray &r, Vector &Point, Vector &Normal, int &sphere_index) const
    {
      bool global_inter = false;
      double min_t = 1E99;
      for (int i = 0; i < objects.size(); i++)
      {
        Vector localPoint, localNormal;
        double local_t;
        bool local_inter = objects[i]->intersect(r, localPoint, localNormal, local_t);
        if (local_inter)
        {
          if (local_t < min_t)
          {
            min_t = local_t;
            Point = localPoint;
            Normal = localNormal;
            sphere_index = i;
          }
          global_inter = true;
        }
      }
      return global_inter;
    }
    
    std::vector<const Object*> objects;
    const Sphere *lumiere;
    double intensity;
};
