//use g++ TD6.cpp -fopenmp dans la terminale

#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include <list>
#include <algorithm>
#include "Vector.h"
#include "Ray.h"
#include "Object.h"
#include "Triangle.h"
#include "Plane.h"
#include "Sphere.h"
#include "BBox.h"
#include "BVH.h"
#include "Geometry.h"
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
Vector get_colour(const Ray&, const Scene&, int, bool show_lights=true);
Vector random_gauss_vect(double, double, double);
Vector random_cos(const Vector&);
Vector random_Phong(const Vector&, double);
std::vector<unsigned char> median_filter(const std::vector<unsigned char>&, int, int);

int main()
{
  //taille image
  int W = 800;
  int H = 800;

  //quantite de rayons a lancer
  const int rays = 120;

  //spheres
  Sphere floor(Vector(0., -1000., 0.), 990, Vector(0., 0., 0.8));
  Sphere sphere(Vector(0., 0., 0.), 10, Vector(1., 0., 0.), 0., 0., false, true);
  /*
  Sphere sphere2(Vector(15., 0., 0.), 5, Vector(1., 0., 0.), 0.2, 1000.);
  Sphere sphere3(Vector(-15., 0., 0.), 5, Vector(1., 0., 0.), 0.8, 1000.);
  Sphere ceiling(Vector(0., 1000., 0.), 940, Vector(0.8, 0., 0.));
  Sphere background(Vector(0., 0., -1000.), 940, Vector(0., 1., 0.));
  Sphere leftwall(Vector(-1000., 0., 0.), 960, Vector(1., 1., 0.));
  Sphere rightwall(Vector(1000., 0., 0.), 960, Vector(1., 0., 0.5));
  */
  
  //triangles

  /*
  Triangle floor(Vector(-1000., -500., -60.), Vector(1000., -500., -60.),
                        Vector(0., -500., 1000.), Vector(0., 0., 0.8),
                        0., 1000., false, false, -1);
  */
  Triangle ceiling(Vector(-1000., 60., -65.), Vector(1000., 60., -65.),
                        Vector(0., 60., 1000.), Vector(0.8, 0., 0.));
  Triangle background(Vector(-1000., -60., -60.), Vector(1000., -60., -60.),
                        Vector(0., 1000., -60.), Vector(0., 1., 0.));
  Triangle leftwall(Vector(-40., -60., -1000.), Vector(-40., -60., 1000.),
                        Vector(-40., 1000., 0.), Vector(1., 1., 0.),
                        0., 1000., false, false, -1);
  Triangle rightwall(Vector(40., -60., -1000.), Vector(40., -60., 1000.),
                        Vector(40., 1000., 0.), Vector(1., 0., 0.5));

  //geometries
  /*
  Geometry girl("Girl/girl.obj", 20.,
                  Vector(0., -9.7, 2.7), Vector(1., 1., 1.),
                  false, false, -1);
  */
  
  //lumiere etendue
  double intensity = 1000000000;
  Sphere lumiere(Vector(0., 40., 20.), 15, Vector(1., 1., 1.), 0, 1000., false, false, true);

  //scene
  Scene scene(&lumiere, intensity);
  scene.addSphere(lumiere);
  scene.addSphere(sphere);
  scene.addSphere(floor);
  /*
  scene.addSphere(sphere2);
  scene.addSphere(sphere3);
  scene.addSphere(ceiling);
  scene.addSphere(background);
  scene.addSphere(leftwall);
  scene.addSphere(rightwall);
  */
  scene.intensity = scene.intensity /
                      (4 * M_PI * scene.lumiere->R * scene.lumiere->R * M_PI);
  scene.addTriangle(ceiling);
  scene.addTriangle(background);
  //scene.addTriangle(floor);
  scene.addTriangle(leftwall);
  scene.addTriangle(rightwall);
  //scene.addGeometry(girl);

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
  std::vector<unsigned char> newimage = median_filter(image, H, W);
  stbi_write_png("imageTD5_fresnel_nomedian.png", W, H, 3, &image[0], 0);
  stbi_write_png("imageTD5_fresnel_median.png", W, H, 3, &newimage[0], 0);

  return 0;
  
}

double Phong_BRDF(const Vector &wi, const Vector &wo, const Vector &N, double exp_phong)
{
  Vector reflected = wo.reflection(N);
  double lobe = std::pow(reflected.dot(wi), exp_phong) * (exp_phong + 2) / (2. * M_PI);

  return lobe;
}

Vector get_colour(const Ray &CurrRay, const Scene &scene, int max_bounces, bool show_lights)
{
  if (max_bounces == 0)
    return Vector(0., 0., 0.);

  //pour obtenir les information de l'intersection avec la sphere
  Vector Normal, Point, colour;
  int index;
  bool inter = scene.intersect(CurrRay, Point, Normal, index, colour);

  Vector sphere_light = Vector(0., 0., 0.);
  if (inter)
  {
    if (scene.objects[index]->is_light)
    {
      sphere_light = show_lights?(scene.intensity * colour):Vector(0., 0., 0.);
    }

    else if (scene.objects[index]->is_mirror)
    {
      Vector reflected = CurrRay.u.reflection(Normal);
      reflected = reflected.getNormalized();
      Ray newray(Point + 0.001 * Normal, reflected);
      sphere_light = get_colour(newray, scene, max_bounces - 1);
    }

    else if (scene.objects[index]->is_transparent)
    {
      double n1 = 1;
      double n2 = 1.3;
      Vector trans_norm = Normal;
      Ray new_trans_ray;
      bool entering = true;
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
        entering = false;
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
        double R0 = sqrt((n1 - n2) / (n1 + n2)),
                R;
        if (entering)
          R = R0 + (1 - R0) * std::pow(1 + CurrRay.u.dot(Normal), 5);
        else
          R = R0 + (1 - R0) * std::pow(1 - trans_vect.dot(Normal), 5);

        if (distribution(generator) < R)
          new_trans_ray = Ray(Point + 0.001 * trans_norm, CurrRay.u.reflection(Normal));
        else
          new_trans_ray = Ray(Point - 0.001 * trans_norm, trans_vect);
      }

      else
        new_trans_ray = Ray(Point + 0.001 * trans_norm, CurrRay.u.reflection(Normal));

      sphere_light = get_colour(new_trans_ray, scene, max_bounces - 1);
    }

    else
    {
      //DEBUT ECLAIRAGE DIRECT
      Vector axe_OP = (Point - scene.lumiere->O).getNormalized();
      Vector random_dir_direct = random_cos(axe_OP);
      Vector random_point = scene.lumiere->R * random_dir_direct + scene.lumiere->O;
      Vector wi = (random_point - Point).getNormalized();
      double dist_light = (random_point - Point).getNorm2();

      /*
        rayon intersection - lumiere. + 0.01 * Normal pour
        eviter les intersections avec le point lui-meme
      */
      Ray shade(Point + 0.001 * Normal, wi);
      int shade_index;
      Vector shade_Normal, shade_Point, shade_colour;

      //pour savoir s'il y a des intersections
      bool shade_inter = scene.intersect(
                                          shade, shade_Point,
                                          shade_Normal, shade_index,
                                          shade_colour
                                        );

      //distances carrees, pas besoin d'obtenir la racine carree
      double dist_point = (shade_Point - Point).getNorm2();

      // 0.99 pour eviter les intersections avec la source de lumiere meme
      if (shade_inter && dist_point < dist_light * 0.99)
      {
        sphere_light = Vector(0., 0., 0.);
      }

      else
      {
        double J = random_dir_direct.dot(-wi) / dist_light,
                prob_func = axe_OP.dot(random_dir_direct) / (M_PI * scene.lumiere->R * scene.lumiere->R);
        Vector BRDF = (
                        (1. - scene.objects[index]->ks) * colour +
                        (
                          Phong_BRDF(wi, CurrRay.u, Normal, scene.objects[index]->phong_exp) *
                          scene.objects[index]->ks
                        ) * colour
                      ) / M_PI;


        sphere_light = scene.intensity * std::max(0., Normal.dot(wi)) * J / prob_func * BRDF;
      }
      //FIN ECLAIRAGE DIRECT

      //DEBUT ECLAIRAGE INDIRECT
      //creation d'une direction aleatoire
      Vector random_dir_indirect,
              reflected_ray = CurrRay.u.reflection(Normal);
      double prob_diffuse = 1 - scene.objects[index]->ks;
      bool diffuse_sample;

      if (distribution(generator) < prob_diffuse)
      {
        diffuse_sample = true;
        random_dir_indirect = random_cos(Normal);
      }
      else
      {
        diffuse_sample = false;
        random_dir_indirect = random_Phong(reflected_ray,
                                            scene.objects[index]->phong_exp);
        if (random_dir_indirect.dot(Normal) < 0 || random_dir_indirect.dot(reflected_ray) < 0)
          return Vector(0., 0., 0.);
      }

      //addition de l'eclairage indirect
      Ray newray(Point + 0.001 * Normal, random_dir_indirect);
      double prob_phong = (scene.objects[index]->phong_exp + 1) / (2 * M_PI) * std::pow(random_dir_indirect.dot(reflected_ray), scene.objects[index]->phong_exp);
      double prob_total = prob_diffuse * Normal.dot(random_dir_indirect) / M_PI +
                          (1. - prob_diffuse) * prob_phong;
      

      if (diffuse_sample)
      {
        sphere_light += get_colour(
                                    newray,
                                    scene, max_bounces - 1,
                                    false
                                  ) * ((Normal.dot(random_dir_indirect) / M_PI / prob_total) * colour);
      }
      else
      {
        sphere_light += get_colour(
                                    newray,
                                    scene, max_bounces - 1,
                                    false
                                  ) * ((
                                        Normal.dot(random_dir_indirect) *
                                        Phong_BRDF(random_dir_indirect, CurrRay.u, Normal, scene.objects[index]->phong_exp) *
                                        scene.objects[index]->ks / prob_total
                                        ) *
                                          colour
                                      );
      }
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
  Vector tangent1_dir;

  double delta_0 = fabs(N[0] - 0.0),
          delta_1 = fabs(N[1] - 0.0),
          delta_2 = fabs(N[2] - 0.0);

  if (delta_0 <= delta_1 && delta_0 <= delta_2)
    tangent1_dir = Vector(0, N[2], -N[1]);
  else if (delta_1 <= delta_0 && delta_1 <= delta_2)
    tangent1_dir = Vector(N[2], 0, -N[0]);
  else
    tangent1_dir = Vector(N[1], -N[0], 0);

  tangent1_dir = tangent1_dir.getNormalized();
  Vector tangent2_dir = tangent1_dir.cross(N);


  return random_local_dir[2] * N +
          random_local_dir[0] * tangent1_dir +
          random_local_dir[1] * tangent2_dir;
}

Vector random_Phong(const Vector &reflected, double exp)
{
  double r1 = distribution(generator),
          r2 = distribution(generator),
          fact = sqrt(1 - std::pow(r2, 2./(exp + 1)));
  Vector random_local_dir = Vector(
                                    cos(2 * M_PI * r1) * fact,
                                    sin(2 * M_PI * r1) * fact,
                                    std::pow(r2, 1./(exp + 1))
                                  );
  Vector tangent1_dir;

  double delta_0 = fabs(reflected[0] - 0.0),
          delta_1 = fabs(reflected[1] - 0.0),
          delta_2 = fabs(reflected[2] - 0.0);

  if (delta_0 <= delta_1 && delta_0 <= delta_2)
    tangent1_dir = Vector(0, reflected[2], -reflected[1]);
  else if (delta_1 <= delta_0 && delta_1 <= delta_2)
    tangent1_dir = Vector(reflected[2], 0, -reflected[0]);
  else
    tangent1_dir = Vector(reflected[1], -reflected[0], 0);

  tangent1_dir = tangent1_dir.getNormalized();
  Vector tangent2_dir = tangent1_dir.cross(reflected);

  return random_local_dir[2] * reflected +
          random_local_dir[0] * tangent1_dir +
          random_local_dir[1] * tangent2_dir;
}

std::vector<unsigned char> median_filter(const std::vector<unsigned char> &im, int H, int W)
{
  std::vector<unsigned char> final_image(W * H * 3, 0);

  for (int i = 0; i < H; i++)
  {
    for (int j = 0; j < W; j++)
    {
      int neighbours[9][2] = {
                                  {i - 1, j - 1},
                                  {i - 1, j},
                                  {i - 1, j + 1},
                                  {i, j - 1},
                                  {i, j},
                                  {i, j + 1},
                                  {i + 1, j - 1},
                                  {i + 1, j},
                                  {i + 1, j + 1}
                              };
      std::vector<int> red_values,
                        green_values,
                        blue_values;

      for (int k = 0; k < 9; k++)
      {
        if (neighbours[k][0] >= 0 && neighbours[k][1] >= 0 &&
            neighbours[k][0] < H && neighbours[k][1] < W)
        {
          red_values.push_back(im[((H - neighbours[k][0] - 1) * W + neighbours[k][1]) * 3 + 0]);
          green_values.push_back(im[((H - neighbours[k][0] - 1) * W + neighbours[k][1]) * 3 + 1]);
          blue_values.push_back(im[((H - neighbours[k][0] - 1) * W + neighbours[k][1]) * 3 + 2]);
        }
      }

      std::sort(red_values.begin(), red_values.end());
      std::sort(green_values.begin(), green_values.end());
      std::sort(blue_values.begin(), blue_values.end());

      final_image[((H - i - 1) * W + j) * 3 + 0] = red_values[ceil(red_values.size() / 2) - 1];
      final_image[((H - i - 1) * W + j) * 3 + 1] = green_values[ceil(green_values.size() / 2) - 1];
      final_image[((H - i - 1) * W + j) * 3 + 2] = blue_values[ceil(blue_values.size() / 2) - 1];
    }
  }

  return final_image;
}

