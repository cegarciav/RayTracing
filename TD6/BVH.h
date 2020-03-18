class BVH
{
  public:
    BVH(){};

    int i_start, i_end;
    BBox bbox;
    BVH *left_child, *right_child;
};
