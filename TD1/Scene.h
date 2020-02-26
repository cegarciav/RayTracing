class Scene
{
  public:
    Scene() {}
    
    void addSphere(const Sphere &s)
    {
    	spheres.push_back(s);
    }

    bool intersect(Ray &r, Vector &Point, Vector &Normal, int &sphere_index)
    {
    	bool global_inter = false;
    	double min_t = 1E99;
    	for (int i = 0; i < spheres.size(); i++)
    	{
    		Vector localPoint, localNormal;
    		double local_t;
		    bool local_inter = spheres[i].intersect(r, localPoint, localNormal, local_t);
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

    Sphere operator[](int index) const
    {
      return spheres[index];
    }

    Sphere operator[](int index)
    {
      return spheres[index];
    }
    
    std::vector<Sphere> spheres;
};
