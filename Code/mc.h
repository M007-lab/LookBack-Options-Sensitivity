#ifndef LOOKBACKOPTION_MC_H
#define LOOKBACKOPTION_MC_H


#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <cmath>
#include "matplotlibcpp.h"
//#include <omp.h>

namespace plt = matplotlibcpp;
struct mean_var {
    mean_var(unsigned n = 0, double sum_x = 0, double sum_xx = 0)
            : sample_size(n), sum_x(sum_x), sum_xx(sum_xx) { }
    double mean() const { return sum_x / (double) sample_size; }
    double var() const { return (sum_xx - sample_size * mean() * mean())
                                / (double) (sample_size-1); }
    double ic_size() const { return 1.96 * std::sqrt(var() / sample_size); }

    mean_var & operator+=(mean_var const & mv) {
        sample_size += mv.sample_size;
        sum_x += mv.sum_x;
        sum_xx += mv.sum_xx;
        return *this;
    }
    friend mean_var operator+(mean_var const & mv1, mean_var const & mv2) {
        return { mv1.sample_size + mv2.sample_size,
                 mv1.sum_x + mv2.sum_x,
                 mv1.sum_xx + mv2.sum_xx };
    }
    friend mean_var operator*(double alpha, mean_var const & mv) {
        return { mv.sample_size, alpha * mv.sum_x, alpha * mv.sum_xx };
    }
    friend std::ostream & operator<<(std::ostream & o, mean_var const & mv) {
        return o << "\tMean: " << mv.mean()
                 << "\t|\tVariance: " << mv.var()
                 << "\t|\tCI: [" << mv.mean()-mv.ic_size() << "," << mv.mean()+mv.ic_size() << "]" ;
    }

protected:
    unsigned sample_size;
    double sum_x, sum_xx;
};


template <typename TDistrib,typename TGen>
mean_var monte_carlo(TDistrib & X,TGen & gen,double discount) {
    double x = X(gen);
    return { 1, discount*x, discount*discount*x*x };
}


template <typename TDistrib, typename TGen>
mean_var monte_carlo(TDistrib & Y, TGen & gen, unsigned sample_size,double discount,std::string title) {
    int n = (int) sample_size/10;
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> ic1(n);
    std::vector<double> ic2(n);
    auto mv = mean_var();
    for (unsigned k = 0; k < sample_size; ++k) {
        mv += monte_carlo(Y,gen,discount);
        if(k%10 ==0){
            int i = k/10;
            x.at(i) =  k;
            y.at(i) = mv.mean();
            ic1.at(i) = mv.mean() - mv.ic_size();
            ic2.at(i) = mv.mean() + mv.ic_size();

        }
    }
    plt::named_plot(title ,x,y,"--");
    // plt::plot(x, ic1,"--"); 
    // plt::plot(x, ic2,"--"); 
    plt::xlabel("Sample Size");
    plt::legend();
    return mv;
}




template <typename TDistrib, typename TGen>
mean_var monte_carlo(TDistrib & Y, TGen & gen, unsigned sample_size,double discount,std::string title,double true_value,double y_min,double y_max) {
    int n = (int) sample_size/10;
    std::vector<double> x(n);
    std::vector<double> y(n);
    std::vector<double> ic1(n);
    std::vector<double> ic2(n);
    std::vector<double> y_star(n);
    std::vector<double> var(n);
    y_star.at(0) = true_value;
    auto mv = mean_var();
    for (unsigned k = 0; k < sample_size; ++k) {
        mv += monte_carlo(Y,gen,discount);
        if(k%10 ==0){
            int i = k/10;
            x.at(i) =  k;
            y.at(i) = mv.mean();
            var.at(i) = mv.var();
            ic1.at(i) = mv.mean() - mv.ic_size();
            ic2.at(i) = mv.mean() + mv.ic_size();
            y_star.at(i) = true_value;

        }
    }
    plt::named_plot(title ,x,y,"--");
    plt::named_plot("True Value",x,y_star);
    // plt::plot(x, ic1, "--"); 
    // plt::plot(x, ic2, "--"); 
    plt::ylim(y_min,y_max);
    plt::xlabel("Sample Size");
    plt::legend();
    return mv;
}



#endif //LOOKBACKOPTION_MC_H
