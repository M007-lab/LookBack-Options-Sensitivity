#ifndef LOOKBACKOPTION_PAYOFFLB_H
#define LOOKBACKOPTION_PAYOFFLB_H

#include "Payoff.h"
#include "Greek.h"

template<typename TRandom>
struct PayoffLb : public Payoff<TRandom>{
    virtual void set_greek(Greek<TRandom>* greek){
        this->greek = greek;
    }

    virtual PayoffLb<TRandom>* antithetic() = 0;
    virtual PayoffLb<TRandom>* shifted() = 0;
    virtual Greek<TRandom>* get_greek(){return this->greek;};
    void set_K(double K){
        this->K = K;
    }

    double get_K() const {
        return K;
    }

    virtual void stop(){
        delete greek;
        greek = nullptr;
    }
    ~PayoffLb(){
        delete greek;
    }

    template<typename TGen>
    double operator()(TGen& gen){
        this->X->template operator()(gen);
        return this->greek->value();
    }

protected:
    Greek<TRandom>* greek;
    double K;
};




#endif 
