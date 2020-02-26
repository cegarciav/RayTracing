//use g++ TD4.cpp -fopenmp dans la terminale

#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

/*
  pour generer des nombres aleatoires
  entre 0 et 1 avec distribution uniforme
*/
static std::default_random_engine generator;
static std::uniform_real_distribution<double> distribution(0.0, 1.0);


/*
  fonction recursive pour obtenir l'intensité
  de la lumière dans un point donné
*/
Vector get_colour(const Ray&, const Scene&, int);
Vector random_gauss_vect(double, double, double);
Vector random_cos(const Vector&);

int main()
{
  //taille image
  int W = 800;
  int H = 800;

  //quantite de rayons a lancer
  const int rays = 500;

  //spheres
  Sphere sphere(Vector(0., 0., 0.), 10, Vector(1., 0., 0.));
  Sphere sphere2(Vector(15., 15., 10.), 10, Vector(0., 1., 0.), false, true);
  Sphere sphere3(Vector(-15., 15., -20.), 10, Vector(0., 0., 1.), true);
  Sphere ceiling(Vector(0., 1000., 0.), 940, Vector(0.8, 0., 0.));
  Sphere background(Vector(0., 0., -1000.), 940, Vector(0., 1., 0.));
  Sphere floor(Vector(0., -1000., 0.), 990, Vector(0., 0., 0.8));
  Sphere leftwall(Vector(-1000., 0., 0.), 960, Vector(1., 1., 0.));
  Sphere rightwall(Vector(1000., 0., 0.), 960, Vector(1., 0., 0.5));

  //lumiere ponctuelle
  double intensity = 1000000000;
  
  //lumiere etendue
  Sphere lumiere(Vector(0., 60., 0.), 15, Vector(1., 1., 1.));

  //scene
  Scene scene(lumiere, intensity);
  scene.addSphere(lumiere);
  scene.addSphere(sphere);
  scene.addSphere(sphere2);
  scene.addSphere(sphere3);
  scene.addSphere(ceiling);
  scene.addSphere(background);
  scene.addSphere(floor);
  scene.addSphere(leftwall);
  scene.addSphere(rightwall);


  //parametres camera
  Vector C(0., 0., 55.);
  double d = W / (2 * sqrt(3) / 3);
  double focus_distance = 55;
  
  std::vector<unsigned char> image(W * H * 3, 0);

  #pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < H; i++)
  {
    for (int j = 0; j < W; j++)
    {
      Vector final_colour(0., 0., 0.);
      for (int ray = 0; ray < rays; ++ray)
      {
        //rayon camera-pixel[i,j]
        double dx_diaphragm = (distribution(generator) - 0.5) * 5.;
        double dy_diaphragm = (distribution(generator) - 0.5) * 5.;

        Vector Xij = random_gauss_vect(j - W / 2, i - H / 2, -d);
        Vector u = Xij.getNormalized();

        Vector destination = C + focus_distance * u;
        Vector origin = C  + Vector(dx_diaphragm, dy_diaphragm, 0);

        Ray Rayij(origin, (destination - origin).getNormalized());
        int index;
        final_colour += get_colour(Rayij, scene, 5) / rays;
      }

      //correction gamma ajoutee      
      image[((H - i - 1) * W + j) * 3 + 0] = 
              std::min(255., std::max(0., std::pow(final_colour[0], 1. / 2.2)));
      image[((H - i - 1) * W + j) * 3 + 1] =
              std::min(255., std::max(0., std::pow(final_colour[1], 1. / 2.2)));
      image[((H - i - 1) * W + j) * 3 + 2] =
              std::min(255., std::max(0., std::pow(final_colour[2], 1. / 2.2)));

    }
  }
  stbi_write_png("imageTD4.png", W, H, 3, &image[0], 0);

  return 0;
}


Vector get_colour(const Ray &CurrRay, const Scene &scene, int max_bounces)
{
  if (max_bounces == 0)
    return Vector(0., 0., 0.);

  //pour obtenir les information de l'intersection avec la sphere
  Vector Normal, Point;
  int index;
  bool inter = scene.intersect(CurrRay, Point, Normal, index);

  Vector sphere_light = Vector(0., 0., 0.);
  if (inter)
  {
    if (scene[index].is_mirror)
    {
      Vector reflected = CurrRay.u - 2 * CurrRay.u.dot(Normal) * Normal;
      reflected = reflected.getNormalized();
      sphere_light = get_colour(
                                  Ray(Point + 0.001 * Normal, reflected),
                                  scene, max_bounces - 1
                                );
    }

    else if (scene[index].is_transparent)
    {
      double n1 = 1;
      double n2 = 1.3;
      Vector trans_norm = Normal;
      /*
        si le produit scalaire entre le rayon et la normal est
        plus grand que zero, alors on est en train de sortir de
        la sphere. Cela veut dire qu'il faut echanger n1 et n2
        et changer la direction du vector normal
      */
      if (CurrRay.u.dot(Normal) > 0)
      {
        double n_aux = n1;
        n1 = n2;
        n2 = n_aux;
        trans_norm = -Normal;
      }
      double normal_coeff = 1. - std::pow((n1 / n2), 2) *
                            (1. - std::pow(CurrRay.u.dot(trans_norm), 2));
      if (normal_coeff > 0)
      {
        Vector trans_vect_t = (n1 / n2) * (
                                            CurrRay.u - 
                                            CurrRay.u.dot(trans_norm) * trans_norm
                                          );
        Vector trans_vect_n = -sqrt(normal_coeff) * trans_norm;
        Vector trans_vect = (trans_vect_t + trans_vect_n).getNormalized();
        sphere_light = get_colour(
                                    Ray(Point - 0.001 * trans_norm, trans_vect),
                                    scene, max_bounces - 1
                                  );
      }
    }

    else
    {
      //DEBUT ECLAIRAGE DIRECT
      Vector axe_OP = (Point - scene.lumiere.O).getNormalized();
      Vector random_dir_direct = random_cos(axe_OP);
      Vector random_point = scene.lumiere.R * random_dir_direct + scene.lumiere.O;
      Vector wi = (random_point - Point).getNormalized();
      double dist_light = (random_point - Point).getNorm2();

      /*
        rayon intersection - lumiere. + 0.01 * Normal pour
        eviter les intersections avec le point lui-meme
      */
      Ray shade(Point + 0.001 * Normal, wi);
      int shade_index;
      Vector shade_Normal, shade_Point;

      //pour savoir s'il y a des intersections
      bool shade_inter = scene.intersect(
                                          shade, shade_Point,
                                          shade_Normal, shade_index
                                        );

      //distances carrees, pas besoin d'obtenir la racine carree
      double dist_point = (shade_Point - Point).getNorm2();

      // 0.99 pour eviter les intersections avec la source de lumiere meme
      if (shade_inter && dist_point < dist_light * 0.99)
      {
        sphere_light = Vector(0., 0., 0.);
      }

      else
        sphere_light = (scene.intensity / (4 * M_PI * dist_light)) *
                        (std::max(0., Normal.dot(wi)) * random_dir_direct.dot(-wi) /
                                      axe_OP.dot(random_dir_direct)) * scene[index].albedo;
      //FIN ECLAIRAGE DIRECT

      //DEBUT ECLAIRAGE INDIRECT
      //creation d'une direction aleatoire
      Vector random_dir_indirect = random_cos(Normal);

      //addition de l'eclairage indirect
      sphere_light += get_colour(
                                  Ray(Point + 0.001 * Normal, random_dir_indirect),
                                  scene, max_bounces - 1
                                ) * scene[index].albedo;
      //FIN ECLAIRAGE INDIRECT

    }

  }

  return sphere_light;
}



Vector random_gauss_vect(double ox, double oy, double z)
{
  double r1 = distribution(generator),
          r2 = distribution(generator),
          R = sqrt(-2 * log(r1)),
          sigma = 0.5;
  //deux nombres generes de maniere aleatoire gaussienne par Box-Muller
  double dx =  R * cos(2 * M_PI * r2) * sigma,
          dy =  R * sin(2 * M_PI * r2) * sigma;

  //0.5 pour se mettre dans le centre du pixel a la place du coin haut-gauche
  return Vector(ox + dx + 0.5, oy + dy - 0.5, z);
}

Vector random_cos(const Vector &N)
{
  double r1 = distribution(generator),
          r2 = distribution(generator);
  Vector random_local_dir = Vector(
                                    cos(2 * M_PI * r1) * sqrt(1 - r2),
                                    sin(2 * M_PI * r1) * sqrt(1 - r2),
                                    sqrt(r2)
                                  );
  Vector tangent1_dir = N.cross(Vector(
                                        distribution(generator) - 0.5,
                                        distribution(generator) - 0.5,
                                        distribution(generator) - 0.5
                                      )
                                );
  tangent1_dir = tangent1_dir.getNormalized();
  Vector tangent2_dir = tangent1_dir.cross(N);

  return random_local_dir[2] * N +
          random_local_dir[0] * tangent1_dir +
          random_local_dir[1] * tangent2_dir;
}
