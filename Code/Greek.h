#ifndef LOOKBACKOPTION_GREEK_H
#define LOOKBACKOPTION_GREEK_H

#include "Payoff.h"

template<typename TRandom>
struct Greek {
    Greek(Payoff<TRandom>* P){
        this->P = P;
    }
    virtual double value() = 0;
    virtual Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) = 0;
protected:
    Payoff<TRandom>* P;
};


#endif //LOOKBACKOPTION_GREEK_H
