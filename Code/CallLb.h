#ifndef LOOKBACKOPTION_LBCALL_H
#define LOOKBACKOPTION_LBCALL_H

#include "PayoffLb.h"
#include "Greek.h"
#include "quantile.h"



// LB call with fixed strike
template<typename TRandom>
struct CallLb : public PayoffLb<TRandom> {
    CallLb(IEngine<TRandom> &engine,double K){
        this->K = K;
        this->X = &engine;
    };

    ~CallLb() {}

    PayoffLb<TRandom>* antithetic() override {
        auto X_anti = this->X->antithetic();
        PayoffLb<TRandom>* payoff_anti = new CallLb<TRandom>(*X_anti,this->K);
        Greek<TRandom>* greek_anti = this->greek->antithetic(payoff_anti);
        payoff_anti->set_greek(greek_anti);
        return payoff_anti;
    }

    PayoffLb<TRandom>* shifted() override {
        IEngine<TRandom>* X_shift = this->X->shifted();
        PayoffLb<TRandom>* payoff_shift= new CallLb<TRandom>(*X_shift,this->K);
        Greek<TRandom>* greek_anti = this->greek->antithetic(payoff_shift);
        payoff_shift->set_greek(greek_anti);
        return payoff_shift;
    }

    ObserverLb<TRandom>* transform(IEngine<TRandom>* X) override {
        return new LBMaximum<TRandom>(*X);
    }

    double price(double x, double r, double sigma, double T) override {
        // BS_Model
        double d = (std::log(x/this->K)+(r+sigma*sigma/2)*T)/(sigma*std::sqrt(T));
        double a1 = normal_cdf(d);
        double a2 = normal_cdf(d-sigma*std::sqrt(T));
        double a3 = normal_cdf(d-2*r*std::sqrt(T)/sigma);
        return x*a1 - std::exp(-r*T)*this->K*a2 + std::exp(-r*T)*sigma*sigma/(2*r)*x*(std::exp(r*T)*a1 -std::pow(x/this->K,-2*r/(sigma*sigma))*a3);
    }

    double delta(double x, double r, double sigma, double T) override {
        double eps = x/100;
        return (price(x+eps,r,sigma,T)-price(x-eps,r,sigma,T))/(2*eps);
    }

    double gamma(double x, double r, double sigma, double T) override {
        double eps = x/100;
        return (price(x+eps,r,sigma,T)+price(x-eps,r,sigma,T) - 2*price(x,r,sigma,T) )/(eps*eps);
    }

    double payoff(double x) override{
        return std::max(x - this->K,0.0);
    }

    double payoff_bin(double x) override{
        return (x > this->K);
    }
};

#endif //LOOKBACKOPTION_LBCALL_H
