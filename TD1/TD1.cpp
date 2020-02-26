#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <math.h>
#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"



int main()
{
  //taille image
  int W = 1024;
  int H = 1024;

  //spheres
  Sphere sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1));
  Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.8, 0, 0));
  Sphere background(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
  Sphere floor(Vector(0, -1000, 0), 990, Vector(0, 0, 0.8));
  Sphere leftwall(Vector(-1000, 0, 0), 960, Vector(1, 1, 0));
  Sphere rightwall(Vector(1000, 0, 0), 960, Vector(1, 0, 0.5));

  //scene
  Scene scene;
  scene.addSphere(sphere);
  scene.addSphere(ceiling);
  scene.addSphere(background);
  scene.addSphere(floor);
  scene.addSphere(leftwall);
  scene.addSphere(rightwall);

  //parametres camera
  Vector C(0, 0, 55);
  int d = W / (2 * sqrt(3) / 3);

  //lumiere
  Vector lumiere(-10, 20, -20);
  double intensity = 300000;
  
  std::vector<unsigned char> image(W * H * 3, 0);
  for (int i = 0; i < H; i++)
  {
    for (int j = 0; j < W; j++)
    {
      //rayon camera-pixel[i,j]
      Vector Xij(j - W / 2, i - H / 2, -d);
      Vector u = Xij.getNormalized();
      Ray Rayij(C, u);

      //pour obtenir les information de l'intersection avec la sphere
      Vector Normal, Point;
      int index;
      bool inter = scene.intersect(Rayij, Point, Normal, index);

      double sphere_light = 0;
      if (inter)
      {
        sphere_light = intensity *
                        Normal.dot((lumiere - Point).getNormalized()) /
                        (lumiere - Point).getNorm2();

      }

      if (sphere_light < 0)
      {
        sphere_light = 0;
      }
      if (sphere_light > 255)
      {
        sphere_light = 255;
      }
      
      image[((H - i - 1) * W + j) * 3 + 0] = sphere_light * scene[index].getAlbedo(0);
      image[((H - i - 1) * W + j) * 3 + 1] = sphere_light * scene[index].getAlbedo(1);
      image[((H - i - 1) * W + j) * 3 + 2] = sphere_light * scene[index].getAlbedo(2);

    }
  }
  stbi_write_png("imageTD1.png", W, H, 3, &image[0], 0);

  return 0;
}
