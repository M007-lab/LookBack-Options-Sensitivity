#ifndef LOOKBACKOPTION_PRICINGENGINE_H
#define LOOKBACKOPTION_PRICINGENGINE_H

#include "IEngine.h"
#include <cmath>

template<typename TRandom>
struct Engine : public IEngine<TRandom>{
    Engine(double x0,double h,double n,TRandom G) : IEngine<TRandom>(x0,h,n,G) {this->path[0] = x0; this->engine = nullptr;};
    Engine(const IEngine<TRandom> &engine) : IEngine<TRandom>(engine) {};
    Engine(IEngine<TRandom> &engine,bool antithetic): IEngine<TRandom>(engine,antithetic) {} ;
    Engine(IEngine<TRandom> &engine,bool antithetic,bool shifted): IEngine<TRandom>(engine,antithetic,shifted) {} ;
    Engine(double x0,IEngine<TRandom> &engine) : IEngine<TRandom>(x0,engine){};

    virtual double get_drift()=0;
    virtual double get_sigma()=0;
    virtual double next(double z,double state);

};

template<typename TRandom>
double Engine<TRandom>::next(double z, double Xk) {
    return Xk + this->get_drift()* Xk *this->h + std::sqrt(this->h)* this->get_sigma()*Xk*z;
}

template<typename TRandom>
struct PricingEngine : public Engine<TRandom>{
    PricingEngine(double x0, double r, double sigma,double h,TRandom G, double n) :
            Engine<TRandom>(x0,h,n,G),r(r),sigma(sigma){};

    PricingEngine(const PricingEngine &engine) : Engine<TRandom>(engine){
        r = engine.get_r();
        sigma = engine.get_sigma();
        this->path[0] = this->init_state;
    }

    PricingEngine(PricingEngine &engine,bool antithetic) : Engine<TRandom>(engine,antithetic){
        r = engine.get_r();
        sigma = engine.get_sigma();
        this->path[0] = this->init_state;
    }

    PricingEngine(PricingEngine &engine,bool antithetic,bool isShifted) : Engine<TRandom>(engine,antithetic,isShifted){
        r = engine.get_r();
        sigma = engine.get_sigma();
        this->path[0] = this->init_state;
    }

    PricingEngine(double x0, PricingEngine &engine) : Engine<TRandom>(x0,engine){
        r = engine.get_r();
        sigma = engine.get_sigma();
        this->path[0] = this->init_state;
    }

    IEngine<TRandom>* copy(double x0) override{
        return new PricingEngine<TRandom>(x0,*this);
    }

    virtual IEngine<TRandom>* antithetic() override{
        return new PricingEngine<TRandom>(*this,1);
    }

    virtual IEngine<TRandom>* shifted() override{
        return new PricingEngine<TRandom>(*this,1,1);
    }

    double get_r() const {return r;}
    double get_sigma() const {return sigma;}

    double next(double z,double state) override;

    double get_drift() override;
    double get_sigma() override;
  


protected:
    double sigma;
    double r;
};


template<typename TRandom>
double PricingEngine<TRandom>::next(double z,double Xk){
    return Xk*std::exp((r-sigma*sigma/2)*this->h + sigma*sqrt(this->h)*z);

}

template<typename TRandom>
double PricingEngine<TRandom>::get_drift(){return r;}
template<typename TRandom>
double PricingEngine<TRandom>::get_sigma(){return sigma;}

#endif 
