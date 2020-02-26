#include <math.h>
#include <iostream>
#include <random>
#include "eigen-master/Eigen/Eigen"
using namespace std;


double p(double);
double p3(double, double, double);
double f(double);
double f3(double, double, double);
double integral(int);
double integral3(int);

class Mvn
{
public:
  Mvn(const Eigen::VectorXd& mu,
      const Eigen::MatrixXd& s)
  {
    mean = mu;
    sigma = s;
  }
  //~Mvn();
  double pdf(const Eigen::VectorXd& x) const
  {
    double n = x.rows();
    double sqrt2pi = std::sqrt(2 * M_PI);
    double quadform  = (x - mean).transpose() * sigma.inverse() * (x - mean);
    double norm = std::pow(sqrt2pi, - n) *
                  std::pow(sigma.determinant(), - 0.5);

    return norm * exp(-0.5 * quadform);
  }
  //Eigen::VectorXd sample(unsigned int nr_iterations = 20) const;
  Eigen::VectorXd mean;
  Eigen::MatrixXd sigma;
};


int main()
{
  cout << integral3(1000000);
  return 0;
}


double p(double x)
{
  double value, exp_val, den_val, sigma = 0.25;
  exp_val = exp(- pow(x, 2) / (2 * pow(sigma, 2)));
  den_val = sigma * sqrt(2 * M_PI);
  value = exp_val / den_val;

  return value;
}

double f(double x)
{
  double value = pow(cos(x), 30);
  
  return value;
}

double integral(int N)
{
  double value = 0, sigma = 0.25;
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0, 1.0);

  for (int i = 0; i <= N; i++)
  {
    double r1 = distribution(generator), r2 = distribution(generator);
    double x1 = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * sigma;
    double current = f(x1) / p(x1);
    value += current;
  }

  return value / N;
}

double f3(double x, double y, double z)
{
  return pow(cos(x * y * z), 30);
}

double p3(double x, double y, double z)
{
  Eigen::MatrixXd sigma(3, 3);
  sigma << 1, 0.45, 0.45,
          0.45, 1, 0.45,
          0.45, 0.45, 1;
  Eigen::VectorXd mean(3);
    mean << 0, 0, 0;
  Eigen::VectorXd vector(3);
    vector << x, y, z;
  Mvn mvn(mean, sigma);
  return mvn.pdf(vector);
}

double integral3(int N)
{
  double value = 0, sigma = 0.25;
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.0, 1.0);

  for (int i = 0; i <= N; i++)
  {
    double r1 = distribution(generator), r2 = distribution(generator);
    double x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * sigma;
    r1 = distribution(generator), r2 = distribution(generator);
    double y = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * sigma;
    r1 = distribution(generator), r2 = distribution(generator);
    double z = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * sigma;
    value += f3(x, y, z) / p3(x, y, z);
  }

  return value / N;
}

