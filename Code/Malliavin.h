#ifndef LOOKBACKOPTION_MALLIAVIN_H
#define LOOKBACKOPTION_MALLIAVIN_H

#include "IEngine.h"
#include "IObserver.h"
#include "PricingEngine.h"
#include "Greek.h"



template<typename TRandom>
struct ExtremeProcess : public IObserver{

    ExtremeProcess(double a,IEngine<TRandom> *engine) : IObserver(*engine) {
        this->engine = engine;
        this->engine->attach(this);
        max_ = 0;
        min_ = 0;
        tau_m = path*0;
        tau_M = path*0;
        this->a = a;
    }
    ~ExtremeProcess(){engine->detach(this);}
    double delta_weight();
    double gamma_weight();
    void update(double z,double u) override;
    void reinit() override;
protected:
    arma::vec tau_m;
    arma::vec tau_M;
    arma::vec Z;
    double DtYs(int k,int l);
    double psiY(int k);
    double psi_primeY(int k);
    IEngine<TRandom>* engine;
    double max_;
    double min_;
    double time_intg; // intégrale de psi(Yt) par rapport à dt
    double sto_intg; // intégrale de psi(Yt) par rapport à dt
    double gamma_time_intg;  // intégrale de psi(Yt)^2 par rapport à dt
    double process; //mu*t + sigma*W_t
    double a;
};

template<typename TRandom>
void ExtremeProcess<TRandom>::update(double z, double u) {
    i++;
    process += engine->get_drift()*engine->get_h() + engine->get_sigma()*std::sqrt(engine->get_h())*z;
    if(process > max_){
        max_ = process;tau_M[i] = i;
    }
    else{tau_M[i] = tau_M[i-1];}

    if(min_ > process){
        min_ = process; tau_m[i] = i;
    }
    else{tau_m[i] = tau_m[i-1];}
    path[i] = max_ - min_;
    time_intg += engine->get_h()* psiY(i-1);
    gamma_time_intg += engine->get_h()* psiY(i-1)*psiY(i-1);
    sto_intg += std::sqrt(engine->get_h())* psiY(i-1)*z;
    Z[i] = z;
}

template<typename TRandom>
double ExtremeProcess<TRandom>::DtYs(int k, int l){
    int t_M = (k <= tau_M[l]);
    int t_m = (k <= tau_m[l]);
    return (t_M - t_m)*engine->get_sigma();
}

template<typename TRandom>
double ExtremeProcess<TRandom>::psiY(int k){
    auto phi = [](double x) {
        if(x < -1){return 0.0;}
        if((x >= -1) and (x < 0)){return std::exp(x*x/(x*x-1));}
        else{return 1.0;}};
    return phi((a - 2*path[k])/a);
}


template<typename TRandom>
double ExtremeProcess<TRandom>::psi_primeY(int k){
    auto phi_prime = [](double x) {
        if((x > -1) and (x < 0)){
            return -2*x*std::exp(x*x/(x*x-1))/((x*x-1)*(x*x-1));
        }
        else{return 0.0;}};
    return -2/a*phi_prime((a - 2*path[k])/a);
}

template<typename TRandom>
double ExtremeProcess<TRandom>::delta_weight() {
    double double_intg = 0;
    double n = path.size()-1;
    double h = engine->get_h();
    for(int k=0;k<n-1;k++){
        double temp = 0;
        for(int l=0;l<k;l++){
            temp+= DtYs(l,k)* psiY(l)*h;
        }
        double_intg += temp * h * psi_primeY(k);
    }
    return (1/engine->get_init_state()) * (1/engine->get_sigma()) *(sto_intg/time_intg + double_intg/(time_intg*time_intg));
}

template<typename TRandom>
double ExtremeProcess<TRandom>::gamma_weight(){
    /*
    double double_intg = 0;
    double n = path.size()-1;
    for(int k=0;k<n-1;k++){
        double temp = 0;
        for(int l=0;l<k;l++){
            temp+= DtYs(l,k)*psiY(l)*engine->get_h();
        }
        double_intg += temp*engine->get_h()*psi_primeY(k)/(engine->get_stochastic_prime(0));
    }*/
    double delta = delta_weight();
    return delta*delta - delta/(engine->get_init_state()) ;
}

template<typename TRandom>
void ExtremeProcess<TRandom>::reinit() {
    i = 0;
    path = 0*path;
    tau_m = path;
    tau_M = path;
    Z = path ;
    max_ = 0;
    min_ = 0;
    process = 0;
    time_intg = 0; // intégrale de psi(Yt) par rapport à dt
    sto_intg = 0;// intégrale de psi(Yt) par rapport à dt
    gamma_time_intg = 0; ;// intégrale de psi(Yt)^2 par rapport à dt
}


template<typename TRandom>
struct Delta_Malv : public Greek<TRandom>{
    Delta_Malv(Payoff<TRandom>* P,double a) : Greek<TRandom>(P){
        this->a = a;
        double x0 = this->P->get_X()->get_init_state();
        this->M = this->P->transform(this->P->get_X());
        this->Y = new ExtremeProcess<TRandom>(a,this->P->get_X());
        this->n = this->P->get_X()->get_n();
    }

    ~Delta_Malv(){
        delete this->Y;
        delete this->M;
    }

    virtual double value() override {
        return (this->P->payoff((*M)[n]))*Y->delta_weight();
    }

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Delta_Malv<TRandom>(P_anti,a);
    }

protected:
    ObserverLb<TRandom>* M;
    double a;
    ExtremeProcess<TRandom>* Y;
    double n;
};

template<typename TRandom>
struct Gamma_Malv : public Delta_Malv<TRandom>{
    Gamma_Malv(Payoff<TRandom>* P,double a) : Delta_Malv<TRandom>(P,a){}


    virtual double value() override {
        return (this->P->payoff(this->M->operator[](this->n)))*this->Y->gamma_weight();
    }

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Gamma_Malv<TRandom>(P_anti,this->a);
    }
};

#endif //LOOKBACKOPTION_MALLIAVIN_H
