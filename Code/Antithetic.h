#ifndef LOOKBACKOPTION_ANTITHETIC_H
#define LOOKBACKOPTION_ANTITHETIC_H

#include "PayoffLb.h"
#include "Martingale.h"

template<typename TRandom>
struct Decorator{
    Decorator(PayoffLb<TRandom>* P) : P(P) {}
    ~Decorator(){};
protected:
    PayoffLb<TRandom>* P;
};

template<typename TRandom>
struct Antithetic : public Decorator<TRandom>{
    Antithetic(PayoffLb<TRandom>* P) : Decorator<TRandom>(P) {
        this->P_anti = this->P->antithetic();
    }

    ~Antithetic(){
        delete P_anti;
    }

    template<typename TGen>
    double operator()(TGen& gen){
        this->P->get_X()->template operator()(gen);
        return 0.5*(this->P->get_greek()->value() + this->P_anti->get_greek()->value());
    }

    double operator()(){
        return 0.5*(this->P->get_greek()->value() + this->P_anti->get_greek()->value());
    }

protected:
    PayoffLb<TRandom>* P_anti;
};

#endif //LOOKBACKOPTION_ANTITHETIC_H
