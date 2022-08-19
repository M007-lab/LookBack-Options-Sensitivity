#ifndef LOOKBACKOPTION_LBOBSERVER_H
#define LOOKBACKOPTION_LBOBSERVER_H

#include "IObserver.h"
#include "IEngine.h"

template<typename TRandom>
struct ObserverLb : public IObserver{
    ObserverLb() : IObserver() {};
    ObserverLb(IEngine<TRandom> &engine) : IObserver(engine) {
        this->engine = &engine;
        this->engine->attach(this);
        prev_state = this->engine->get_init_state();
    }
    ~ObserverLb(){engine->detach(this);}
    virtual void update(double z,double u) = 0;
    virtual void reinit();

protected:
    IEngine<TRandom>* engine;
    double prev_state;
    double state;
};

template<typename TRandom>
void ObserverLb<TRandom>::reinit() {
    prev_state = engine->get_init_state();
    path = 0*path;
    i = 0;
    path[0] = prev_state;
    state = prev_state;
}

template<typename TRandom>
struct LBMaximum : public ObserverLb<TRandom>{
    LBMaximum() : ObserverLb<TRandom>() {};
    LBMaximum(IEngine<TRandom> &engine) : ObserverLb<TRandom>(engine){}
    void update(double z,double u) override;
};

template<typename TRandom>
void LBMaximum<TRandom>::update(double z, double u) {
    this->i++;
    double y = this->engine->get_state();
    double sig2 = this->engine->get_sigma() * this->prev_state;
    sig2 *= sig2;
    double max_ = 0.5*(this->prev_state + y + std::sqrt((this->prev_state - y)*(this->prev_state - y) - 2*this->engine->get_h()*sig2*std::log(1-u)));
    this->state = std::max(max_,this->state);
    this->path[this->i] = this->state;
}

template<typename TRandom>
struct LBMinimum : public ObserverLb<TRandom>{
    LBMinimum(IEngine<TRandom> &engine) : ObserverLb<TRandom>(engine){};
    void update(double z,double u) override;
};

template<typename TRandom>
void LBMinimum<TRandom>::update(double z, double u) {
    this->i++;
    double y = this->engine->get_state();
    double sig2 = this->engine->get_sigma() * this->prev_state ;
    sig2 *= sig2;
    double min_ = 0.5*(this->prev_state + y - std::sqrt((this->prev_state - y)*(this->prev_state - y) - 2*this->engine->get_h()*sig2*std::log(1-u)));
    this->state = std::max(min_,this->state);
    this->path[this->i] = this->state;
}
#endif //LOOKBACKOPTION_LBOBSERVER_H
