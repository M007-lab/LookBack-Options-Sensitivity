#ifndef LOOKBACKOPTION_PAYOFF_H
#define LOOKBACKOPTION_PAYOFF_H

#include "IEngine.h"
#include "ObserverLb.h"

template<typename TRandom>
struct Payoff{
    virtual ~Payoff()=0;
    virtual ObserverLb<TRandom>* transform(IEngine<TRandom>* X) = 0;
    virtual double payoff(double x)=0;
    virtual double payoff_bin(double x)=0;
    virtual double price(double x, double r, double sigma, double T) = 0;
    virtual double delta(double x, double r, double sigma, double T) = 0;
    virtual double gamma(double x, double r, double sigma, double T) = 0;
    IEngine<TRandom>* get_X() {return this->X;};
protected:
    IEngine<TRandom>* X;
};

template<typename TRandom>
Payoff<TRandom>::~Payoff(){};

#endif //LOOKBACKOPTION_PAYOFF_H
