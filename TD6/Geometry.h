class Geometry : public Object
{
  public:
    Geometry(
              const char *obj, double scaling,
              const Vector &offset, const Vector &couleur,
              double ks_value=0., double phong_exponent=1000.,
              bool mirror=false, bool transparent=false,
              int new_normal_sign=1
            )
    {
      albedo = couleur;
      ks = ks_value;
      phong_exp = phong_exponent;
      is_mirror = mirror;
      is_transparent = transparent;
      normal_sign = new_normal_sign;
      char matfile[255];

      FILE *f;
      f = fopen(obj, "r");
      int curGroup = -1;
      char line[255];
      while (fgets(line, 255, f))
      {
        if (line[0] == 'm' && line[1] == 't' && line[2] == 'l')
          sscanf(line, "mtllib %s\n", matfile);

        if (line[0] == 'u' && line[1] == 's')
          curGroup++;

        if (line[0] == 'v' && line[1] == ' ')
        {
          double i0, i1, i2;
          sscanf(line, "v %lf %lf %lf\n", &i0, &i2, &i1);
          Vector vec = scaling * Vector(i0, i1, i2) + offset;
          vertices.push_back(vec);
        }

        if (line[0] == 'v' && line[1] == 'n')
        {
          double i0, i1, i2;
          sscanf(line, "vn %lf %lf %lf\n", &i0, &i2, &i1);
          Vector vec(i0, i1, i2);
          normals.push_back(vec);
        }

        if (line[0] == 'v' && line[1] == 't')
        {
          double i0, i1;
          sscanf(line, "vt %lf %lf\n", &i0, &i1);
          Vector vec(i0, i1);
          uvs.push_back(vec);
        }

        if (line[0] == 'f' && line[1] == ' ')
        {
          int i0, i1, i2,
              j0, j1, j2,
              k0, k1, k2;
          faceGroup.push_back(curGroup);
          int success_n = sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                                    &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2);
          if (success_n == 9)
          {
            faces.push_back(i0 - 1);
            faces.push_back(i1 - 1);
            faces.push_back(i2 - 1);

            uvIDs.push_back(j0 - 1);
            uvIDs.push_back(j1 - 1);
            uvIDs.push_back(j2 - 1);

            normalIDs.push_back(k0 - 1);
            normalIDs.push_back(k1 - 1);
            normalIDs.push_back(k2 - 1);
          }
        }
      }

      fclose(f);


      char folder[255];
      sscanf(obj, "%[^/]s%*s", folder);
      f = fopen((std::string(folder) + std::string("/") + std::string(matfile)).c_str(), "r");

      while (fgets(line, 255, f))
      {
        if (line[0] == 'm' && line[4] == 'K' && line[5] == 'd')
        {
          char texture[255];
          sscanf(line, "map_Kd %s\n", texture);
          add_texture((std::string(folder) + std::string("/") + std::string(texture)).c_str());
        }
      }

      fclose(f);

      build_bvh(&bvh, 0, faces.size() / 3);
    }

    bool intersect(const Ray &r, Vector &Point, Vector &Normal, double &t, Vector &colour) const
    {
      bool has_inter = false;
      t = 1E99;
      if (!bvh.bbox.intersect(r))
      {
        return false;
      }

      std::list<const BVH*> bvh_list;
      bvh_list.push_front(&bvh);

      while(!bvh_list.empty())
      {
        const BVH *current = bvh_list.front();
        bvh_list.pop_front();

        if (current->left_child && current->left_child->bbox.intersect(r))
        {
          bvh_list.push_back(current->left_child);
        }
        if (current->right_child && current->right_child->bbox.intersect(r))
        {
          bvh_list.push_back(current->right_child);
        }

        if (!current->left_child)
        {
          for (int i = current->i_start; i < current->i_end; i++)
          {
            int v0 = faces[i * 3],
                v1 = faces[i * 3 + 1],
                v2 = faces[i * 3 + 2];
            Triangle current_triangle(
                                        vertices[v0], vertices[v1],
                                        vertices[v2], albedo,
                                        is_mirror, is_transparent,
                                        normal_sign
                                      );

            Vector localPoint, localNormal, localColour;
            double local_t, alpha, beta, gamma;
            if (current_triangle.intersect(r, localPoint, localNormal, local_t,
                                              alpha, beta, gamma, localColour))
            {
              has_inter = true;
              if (local_t < t)
              {
                t = local_t;
                Point = localPoint;
                if (normals.size() > 0)
                {
                  Normal = alpha * normals[normalIDs[i * 3]] +
                            beta * normals[normalIDs[i * 3 + 1]] +
                            gamma * normals[normalIDs[i * 3 + 2]];
                  Normal = Normal.getNormalized();
                }
                else
                  Normal = localNormal;

                if (textures.size() > 0)
                {
                  int textureId = faceGroup[i];

                  int x = (alpha * uvs[uvIDs[i * 3]][0] +
                          beta * uvs[uvIDs[i * 3 + 1]][0] +
                          gamma * uvs[uvIDs[i * 3 + 2]][0])
                          * (w[textureId] - 1),

                      y = (alpha * uvs[uvIDs[i * 3]][1] +
                          beta * uvs[uvIDs[i * 3 + 1]][1] +
                          gamma * uvs[uvIDs[i * 3 + 2]][1])
                          * (h[textureId] - 1);

                  double cr = (textures[textureId][(y * w[textureId] + x) * 3]) / 255.,
                          cg = (textures[textureId][(y * w[textureId] + x) * 3 + 1]) / 255.,
                          cb = (textures[textureId][(y * w[textureId] + x) * 3 + 2]) / 255.;


                  colour = Vector(cr, cg, cb);
                }
                else
                  colour = localColour;
              }
            }
          }
        }
      }

      return has_inter;      
    }

    void add_texture(const char *filename)
    {
      textures.resize(textures.size() + 1);
      w.resize(w.size() + 1);
      h.resize(h.size() + 1);

      FILE *f;
      f = fopen(filename, "rb");
      unsigned char info[54];
      fread(info, sizeof(unsigned char), 54, f); //read header 54-byte

      w[w.size() - 1] = *(int*)&info[18];
      h[h.size() - 1] = *(int*)&info[22];

      int size = 3 * w[w.size() - 1] * h[h.size() - 1];
      textures[textures.size() - 1].resize(size);
      fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f);
      fclose(f);

      for (int i = 0; i < size; i += 3)
      {
        std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
      }
    }

    BBox build_bbox(int i0, int i1)
    {
      BBox newbbox;
      double min_vec[3] = {
                            vertices[faces[i0 * 3]][0],
                            vertices[faces[i0 * 3]][1],
                            vertices[faces[i0 * 3]][2]
                          },
              max_vec[3] = {
                              vertices[faces[i0 * 3]][0],
                              vertices[faces[i0 * 3]][1],
                              vertices[faces[i0 * 3]][2]
                            };

      for (int triang = i0; triang < i1; triang++) //indice de triangle
      {
        for (int vert = 0; vert < 3; vert++) //indice de sommet
        {
          for (int coord = 0; coord < 3; coord++) //indice de coordonnee
          {
            min_vec[coord] = std::min(
                                        vertices[faces[triang * 3 + vert]][coord],
                                        min_vec[coord]
                                      );
            max_vec[coord] = std::max(
                                        vertices[faces[triang * 3 + vert]][coord],
                                        max_vec[coord]
                                      );
          }
        }
      }

      newbbox.bmin = Vector(min_vec[0], min_vec[1], min_vec[2]);
      newbbox.bmax = Vector(max_vec[0], max_vec[1], max_vec[2]);

      return newbbox;
    }

    void build_bvh(BVH *node, int i0, int i1)
    {
      node->bbox = build_bbox(i0, i1);
      node->i_start = i0;
      node->i_end = i1;
      node->left_child = NULL;
      node->right_child = NULL;

      Vector diag = node->bbox.bmax - node->bbox.bmin;
      int split_dim;

      if (diag[0] > diag[1])
      {
        if (diag[0] > diag[2])
          split_dim = 0;
        else
          split_dim = 2;
      }
      else
      {
        if (diag[1] > diag[2])
          split_dim = 1;
        else
          split_dim = 2;
      }

      double split_val = node->bbox.bmin[split_dim] + diag[split_dim] * 0.5;
      int pivot = i0 - 1;
      for (int i = i0; i < i1; i++)
      {
        double center_split = (
                                vertices[faces[i * 3]][split_dim] +
                                vertices[faces[i * 3 + 1]][split_dim] +
                                vertices[faces[i * 3 + 2]][split_dim]) / 3;

        if (center_split < split_val)
        {
          pivot++;
          std::swap(faces[i * 3], faces[pivot * 3]);
          std::swap(faces[i * 3 + 1], faces[pivot * 3 + 1]);
          std::swap(faces[i * 3 + 2], faces[pivot * 3 + 2]);

          if (normals.size() > 0)
          {
            std::swap(normalIDs[i * 3], normalIDs[pivot * 3]);
            std::swap(normalIDs[i * 3 + 1], normalIDs[pivot * 3 + 1]);
            std::swap(normalIDs[i * 3 + 2], normalIDs[pivot * 3 + 2]);
          }

          if (uvs.size() > 0)
          {
            std::swap(uvIDs[i * 3], uvIDs[pivot * 3]);
            std::swap(uvIDs[i * 3 + 1], uvIDs[pivot * 3 + 1]);
            std::swap(uvIDs[i * 3 + 2], uvIDs[pivot * 3 + 2]);

            std::swap(faceGroup[i], faceGroup[pivot]);
          }
        }
      }

      if (pivot <= i0 || pivot >= i1)
        return;

      node->left_child = new BVH();
      build_bvh(node->left_child, i0, pivot);
      node->right_child = new BVH();
      build_bvh(node->right_child, pivot, i1);
    }

    //avec les indices des groupes existants
    std::vector<int> faceGroup;
    /*
      indices des sommets des triangles. Les
      trois premiers elements sont les indices
      du premier triangle, etc
    */
    std::vector<int> faces;

    std::vector<int> normalIDs;
    std::vector<int> uvIDs;
    ///liste de sommets (points en trois dimensions)
    std::vector<Vector> vertices;
    //normales par sommet
    std::vector<Vector> normals;
    //vector avec les coordonnees des textures
    std::vector<Vector> uvs;

  private:
    BVH bvh;
    int normal_sign;
    std::vector<std::vector<unsigned char>> textures;
    std::vector<int> w, h;
};
