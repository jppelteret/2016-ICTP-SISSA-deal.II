
#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#  include <sys/time.h>
#  include <sys/resource.h>
#include "/home/kronbichler/Programme/compilers/iaca-lin64/include/iacaMarks.h"

template <typename Number>
Number horner_evaluation (const std::vector<Number> &coefficients,
                          const Number               x)
{
  const unsigned int N = coefficients.size()-1;
  Number res = coefficients[N];
  for (int j=N-1; j>=0; --j)
    res = coefficients[j] + x * res;
  return res;
}


template <typename Number>
void
run_bench(const unsigned int degree)
{
  const unsigned long long int n_operations = 1000000000ull;
  const unsigned long long int repetitions = n_operations / degree;
  std::vector<Number> coefficients(degree+1);
  std::vector<Number> results(128);
  std::vector<Number> input_x(repetitions/128+1);
  for (unsigned int i=0; i<coefficients.size(); ++i)
    coefficients[i] = -1 + 2.*(double)rand()/RAND_MAX;
  for (unsigned int i=0; i<input_x.size(); ++i)
    input_x[i] = (double)rand()/RAND_MAX;

  struct timeval wall_timer;
  gettimeofday(&wall_timer, NULL);
  double start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

  for (unsigned int i=0; i<repetitions; ++i)
    {
      const unsigned long long int xind = i/128;
      const unsigned int rind = i%128;
      results[rind] += horner_evaluation(coefficients, input_x[xind]);
    }

  gettimeofday(&wall_timer, NULL);
  const double compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  Number check = 0;
  for (unsigned int i=0; i<results.size(); ++i)
    check += results[i];

  std::cout << std::setw(10) << degree << " " << std::setw(12) << 2.e-9*(double)repetitions*degree/compute_time << " " << std::setw(10) << check << std::endl;
}


unsigned int memory_function_1(const unsigned int n,
                               const std::vector<unsigned int> &numbers,
                               const unsigned int jump)
{
  double * t = new double[n];

  t[jump] = numbers[0];
  double return_value = t[jump] * 2;

  delete[] t;
  return (unsigned int)return_value;
}


unsigned int memory_function_2(const unsigned int n,
                               const std::vector<unsigned int> &numbers,
                               const unsigned int jump)
{
  return (unsigned int) std::pow((double)numbers[jump]+0.2, 2.123);
}

void run_memory_bench()
{
  const unsigned long long int n_operations = 100000000ull;
  const unsigned int length = 4000;
  std::vector<unsigned int> lengths(length);
  for (unsigned int i=0; i<lengths.size(); ++i)
    lengths[i] = rand()%lengths.size() + 8;

  unsigned int res = 0;

  struct timeval wall_timer;
  gettimeofday(&wall_timer, NULL);
  double start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;
  for (unsigned int i=0; i<n_operations; )
    for (unsigned int j=0; j<length; ++j, ++i)
      res += memory_function_1(lengths[j], lengths, lengths[j] / 4);

  gettimeofday(&wall_timer, NULL);
  const double compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);

  std::cout << "Time per allocation: " << std::setw(12) << compute_time / n_operations << " " << std::setw(10) << res << std::endl;
}


int main (int argc, char**argv)
{
  typedef double Number;

  if (false)
    {
      // read the polynomial degree
      const unsigned long long int n_operations = 1000000000ull;
      const unsigned int degree = atoi(argv[1]);
      const unsigned long long int repetitions = n_operations / degree;

      std::vector<Number> coefficients(degree+1);
      std::vector<Number> results(128);
      std::vector<Number> input_x(repetitions/128);
      for (unsigned int i=0; i<coefficients.size(); ++i)
        coefficients[i] = -1 + 2.*(double)rand()/RAND_MAX;
      for (unsigned int i=0; i<input_x.size(); ++i)
        input_x[i] = (double)rand()/RAND_MAX;

      struct timeval wall_timer;
      gettimeofday(&wall_timer, NULL);
      double start = wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec;

      for (unsigned int i=0; i<repetitions; ++i)
        {
          const unsigned long long int xind = i/128;
          const unsigned int rind = i%128;
          results[rind] += horner_evaluation(coefficients, input_x[xind]);
        }

      gettimeofday(&wall_timer, NULL);
      const double compute_time = (wall_timer.tv_sec + 1.e-6 * wall_timer.tv_usec - start);
      std::cout << "Computation time: " << compute_time << "s" << std::endl
                << "GFLOPs:           " << 2.e-9*(double)n_operations/compute_time << std::endl;

      Number check = 0;
      for (unsigned int i=0; i<results.size(); ++i)
        check += results[i];
      std::cout << "Result check      " << check << std::endl;
    }
  else if (true)
    {
      for (unsigned int i=1; i<4; ++i)
        run_bench<Number>(i);
      unsigned int factor = 1;
      while (factor < 100)
        {
          for (unsigned int i=4; i<8; ++i)
            run_bench<Number>(factor*i);
          factor *= 2;
        }
    }

  run_memory_bench();

}
