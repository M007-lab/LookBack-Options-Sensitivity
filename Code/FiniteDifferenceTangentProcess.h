#ifndef LOOKBACKOPTION_FINITEDIFFERENCEANDTP_H
#define LOOKBACKOPTION_FINITEDIFFERENCEANDTP_H

#include "Greek.h"

template<typename TRandom>
struct Delta_FD : public Greek<TRandom>{
    Delta_FD(Payoff<TRandom>* P) : Greek<TRandom>(P){
        double x0 = this->P->get_X()->get_init_state();
        this->eps = x0/100;
        this->X_e1 = this->P->get_X()->copy(x0+this->eps);
        this->X_e2 = this->P->get_X()->copy(x0-this->eps);
        this->M_e1 = this->P->transform((this->X_e1));
        this->M_e2 = this->P->transform((this->X_e2));
        this->n = this->P->get_X()->get_n();
    }

    virtual Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Delta_FD<TRandom>(P_anti);
    }

    void set_eps(double e) {
        double x0 = this->P->get_X()->get_init_state();
        this->eps = e;
        this->X_e1->set_init_state(x0+e);
        this->X_e2->set_init_state(x0-e);
        this->X_e1->reinit();
        this->X_e2->reinit();
    }

    ~Delta_FD(){
        delete this->M_e2;
        delete this->M_e1;
        delete this->X_e2;
        delete this->X_e1;
    }

    virtual double value() override {
        return (this->P->payoff(M_e1->operator[](n)) - this->P->payoff(M_e2->operator[](n)))/(2*eps);
    }
protected:
    double eps;
    IEngine<TRandom>* X_e1;
    IEngine<TRandom>* X_e2;
    ObserverLb<TRandom>* M_e1;
    ObserverLb<TRandom>* M_e2;
    double n;
};


template<typename TRandom>
struct Delta_DFD : public Delta_FD<TRandom>{
    Delta_DFD(Payoff<TRandom>* P, double c) : Delta_FD<TRandom>(P){
        this->c = c;
        this->set_eps(c*std::pow(i,1/4 + 0.01));
    }

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Delta_DFD<TRandom>(P_anti,c);
    }

    double value() override {
        double result = (this->P->payoff(this->M_e1->operator[](this->n)) - this->P->payoff(this->M_e2->operator[](this->n)))/(2*this->eps);
        i++;
        this->set_eps(c*std::pow(i,1/4 + 0.01));
        return result;
    }

protected:
    unsigned i = 1;
    double c;
};

template<typename TRandom>
struct Delta_TP : public Delta_FD<TRandom>{
    Delta_TP(Payoff<TRandom>* P) : Delta_FD<TRandom>(P){
        M = this->P->transform((this->P->get_X()));
    }

    ~Delta_TP(){ delete M;}
    double value() override {
        return this->P->payoff_bin(M->operator[](this->n))*(this->M_e1->operator[](this->n) - this->M_e2->operator[](this->n))/(2*this->eps);;
    }

    virtual Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Delta_TP<TRandom>(P_anti);
    }

protected:
    ObserverLb<TRandom>* M;
};

template<typename TRandom>
struct Gamma_FD : public Delta_TP<TRandom>{
    Gamma_FD(Payoff<TRandom>* P) : Delta_TP<TRandom>(P){}

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Gamma_FD<TRandom>(P_anti);
    }

    double value() override {
        return  (this->P->payoff(this->M_e1->operator[](this->n)) + this->P->payoff(this->M_e2->operator[](this->n)) - 2*this->P->payoff(this->M->operator[](this->n)))/(this->eps*this->eps);
    }


};


template<typename TRandom>
struct Gamma_FD_TP : public Delta_FD<TRandom>{
    Gamma_FD_TP(Payoff<TRandom>* P) : Delta_FD<TRandom>(P){
        x0 = this->P->get_X()->get_init_state();
    }

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Gamma_FD_TP<TRandom>(P_anti);
    }

    double value() override {
        return (this->P->payoff_bin(this->M_e1->operator[](this->n))*this->M_e1->operator[](this->n)/(this->M_e1->operator[](0)) - this->P->payoff_bin(this->M_e2->operator[](this->n))*this->M_e2->operator[](this->n)/(this->M_e2->operator[](0)))/(2*this->eps);
    }


protected:
    double x0;
};

#endif //LOOKBACKOPTION_FINITEDIFFERENCEANDTP_H
