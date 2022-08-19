#ifndef LOOKBACKOPTION_IENGINE_H
#define LOOKBACKOPTION_IENGINE_H


#include "IObserver.h"
#include "IObservable.h"

class NotImplementedException
        : public std::exception {

public:

    // Construct with given error message:
    NotImplementedException(const char * error = "Functionality not yet implemented!")
    {
        errorMessage = error;
    }

    // Provided for compatibility with std::exception.
    const char * what() const noexcept
    {
        return errorMessage.c_str();
    }

private:

    std::string errorMessage;
};

template<typename TRandom>
struct IEngine : public IObserver,IObservable{
    IEngine(double x0,double h,double n,TRandom G) : IObserver(n),n_max(n),h(h),G(G),init_state(x0),state(x0){path[0] = x0; engine = nullptr;};
    IEngine(const IEngine &engine);
    IEngine(IEngine &engine,bool antithetic);
    IEngine(IEngine &engine,bool antithetic,bool isShifted);
    IEngine(double x0,IEngine &engine);
    ~IEngine();
    double get_h() const;
    double get_n() const;
    double get_init_state() const;
    double get_theta() const;
    void set_theta(double theta);
    void set_init_state(double x);
    virtual double get_state() const;
    TRandom get_G() const;

    virtual double get_drift() = 0; 
    virtual double get_sigma() = 0;


    virtual double next(double z,double Xk) = 0;
    virtual void update(double z, double u) override;
    virtual void reinit() override;
    void reinitObservers(){
        this->observers = std::list<IObserver*>();
    };

    virtual IEngine<TRandom>* copy(double x0) = 0;

    virtual IEngine<TRandom>* antithetic() = 0;

    virtual IEngine<TRandom>* shifted() = 0;

    template<typename TGen>
    void operator()(TGen & gen);

protected:
    IEngine<TRandom>* engine;
    double theta;
    bool isShifted = 0;
    double init_state;
    double n_max;
    double state;
    double h;
    TRandom G;
    bool isAntithetic = 0;
};

template<typename TRandom>
IEngine<TRandom> ::IEngine(const IEngine<TRandom> &engine): IObserver(engine){
    n_max = engine.get_n();
    h  = engine.get_h();
    G = engine.get_G();
    state = engine.get_state();
    isAntithetic = 0;
    init_state = engine.get_init_state();
    this->engine = nullptr;
    this->observers = std::list<IObserver*>();
}

template<typename TRandom>
IEngine<TRandom>::IEngine(IEngine<TRandom> &engine, bool antithetic) : IEngine(engine) {
    isAntithetic = antithetic;
    this->engine = &engine;
    this->engine->attach(this);
}

template<typename TRandom>
IEngine<TRandom>::IEngine(IEngine<TRandom> &engine, bool antithetic,bool isShifted) : IEngine(engine) {
    isAntithetic = antithetic;
    this->isShifted = isShifted;
    theta = 0;
    this->engine = &engine;
    this->engine->attach(this);
}


template<typename TRandom>
IEngine<TRandom>::IEngine(double x0,IEngine<TRandom> &engine) : IEngine(engine){
    state = x0;
    init_state = x0;
    this->engine = &engine;
    this->engine->attach(this);
}


template<typename TRandom>
IEngine<TRandom>::~IEngine(){
    if(engine != nullptr)
        engine->detach(this);
}

template<typename TRandom>
double IEngine<TRandom>::get_init_state() const {return this->init_state;}

template<typename TRandom>
void IEngine<TRandom>::set_init_state(double x){
    this->init_state = x;
    this->path[0] = x;
}

template<typename TRandom>
double IEngine<TRandom>::get_state() const {
    return this->state;
}

template<typename TRandom>
double IEngine<TRandom>::get_n() const{return n_max;}

template<typename TRandom>
double IEngine<TRandom>::get_h() const{return h;}

template<typename TRandom>
TRandom IEngine<TRandom>::get_G() const{return G;}

template<typename TRandom>
double IEngine<TRandom>::get_theta() const{
    return this->theta;
};

template<typename TRandom>
void IEngine<TRandom>::set_theta(double theta){
    this->theta = theta;
};

/*
template<typename TRandom>
void IEngine<TRandom>::next(double z){
    state += get_drift(state)*h + get_stochastic(state)*std::sqrt(h)*z;
}*/

template<typename TRandom>
void IEngine<TRandom>::reinit(){
    i = 0;
    path = 0*path;
    path[0] = init_state;
    state = init_state;
    auto iterator = observers.begin();
    while (iterator != observers.end()) {
        (*iterator)->reinit();
        ++iterator;}
    ;}


template<typename TRandom>
void IEngine<TRandom>::update(double z, double u){
    i++;
    double shift = (isShifted and i == 1) ? theta : 0;
    double z_ = isAntithetic ? -z : z;
    this->state = next(z_+shift,this->state);
    path[i] = state;
    notify(z_+shift,u);
}

template<typename TRandom>
template<typename TGen>
void IEngine<TRandom>::operator()(TGen & gen){
    this->reinit();
    while(i<n_max){
        double z = G(gen);
        double u = std::uniform_real_distribution<>(0,1)(gen);
        update(z,u);
    }
}



#endif //LOOKBACKOPTION_IENGINE_H
