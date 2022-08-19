#ifndef LOOKBACKOPTION_MARTINGALEWEIGHT_H
#define LOOKBACKOPTION_MARTINGALEWEIGHT_H


template<typename TRandom>
struct Brownian : public IObserver{
    Brownian(IEngine<TRandom> *engine) : IObserver(*engine) {
        path[0] = 0;
        this->engine = engine;
        this->engine->attach(this);
    }


    ~Brownian(){engine->detach(this);}
    void update(double z,double u) override{
        i++;
        path[i] = z;
    };

protected:
    IEngine<TRandom>* engine;
};

template<typename TRandom>
struct TangentProcess : public IObserver{
    TangentProcess(IEngine<TRandom> *engine,bool isBS) : IObserver(*engine) {
        path[0] = 1;
        this->engine = engine;
        this->engine->attach(this);
        this->state = 1;
        isBlackScholes = isBS;
    }

    TangentProcess(PricingEngine<TRandom> &engine) : IObserver(engine) {
        path[0] = 1;
        this->engine = &engine;
        this->engine->attach(this);
        this->state = 1;
        isBlackScholes = 1;
    }
    ~TangentProcess(){engine->detach(this);}
    void update(double z,double u) override;

protected:
    IEngine<TRandom>* engine;
    bool isBlackScholes = 0;
    double state;
};

template<typename TRandom>
void TangentProcess<TRandom>::update(double z, double u) {
    i++;
    double X = engine->get_state();
    if(isBlackScholes){
        path[i] = X/engine->get_init_state();
    }
    else{
        this->state += engine->get_drift()*this->state + engine->get_sigma()*this->state*z;
        path[i] = this->state;
    }
}


template<typename TRandom>
struct Delta_Martingale : public Greek<TRandom>{
    Delta_Martingale(Payoff<TRandom>* P,double a,double lambda,bool isBs) : Greek<TRandom>(P){
        this->a = a;
        x0 = this->P->get_X()->get_init_state();
        this->lambda = lambda;
        this->M = this->P->transform(this->P->get_X());
        this->Y = new TangentProcess<TRandom>(this->P->get_X(),isBs);
        this->dW = new Brownian<TRandom>(this->P->get_X());
        n = this->P->get_X()->get_n();
        H = arma::vec(n/2);
        H[0] = 0;
        this->isBs = isBs;
    }

    ~Delta_Martingale(){
        delete M;
        delete Y;
        delete dW;
    }

    virtual double value() override {
        return (this->P->payoff(M->operator[](n)))*mrtg_weight();
    }

    Greek<TRandom>* antithetic(Payoff<TRandom>* P_anti) override {
        return new Delta_Martingale<TRandom>(P_anti,a,lambda,isBs);
    }

    double adapted_process(int i){
        double h = this->P->get_X()->get_h();
        double T = n*h;
        double x = this->P->get_X()->operator[](i);
        double a_ =  a - a/10000;
        return 1 / (std::pow(std::min(std::abs(x-(x0-a_)),std::abs(x-(x0+a_))),2)*(T/2-i*h));
    }

    int get_tau(){
        H *=0;
        int n_2 = (int) n/2;
        int i = 0, tau = 0;
        double h = this->P->get_X()->get_h();
        double itg = 0;
        while(itg < 1 - 1/100000  and i < n_2){
            // temp = integral;
            // H[i] = adapted_process(i);
            // integral += h*(H[i]+adapted_process(i+1))/2;
           
            itg += (h/lambda) * adapted_process(i);
            i++;
        }
        tau = i - 1 ;
        for (int i=0;i<=n; i++){
            if (i < tau){ H[i] = (1/lambda) * adapted_process(i); }
            else {H[i] = 0.;}
        }
        return tau;
        // if(integral < lambda){
        //     H *= 1/integral;
        // }
        // else if(temp == 0)
        //     H[0] = 1;
        // else
        //     H *= 1/temp;
        // return std::min(i-1,n_2);
    }

    double mrtg_weight(){
        int tau = get_tau();
        // std::cout << "tau = " << tau << "n=" <<n <<"H[n]" << H[n] << std::endl;
        double weight = 0;
        double h = this->P->get_X()->get_h();
        for(int i = 0;i<=n;i++){
            double x = this->P->get_X()->operator[](i);
            double y = this->Y->get_path()[i];
            double sig = this->P->get_X()->get_sigma();
            weight += std::sqrt(h)*this->dW->operator[](i)*H[i]*y / sig;
        }
        
        return weight;
    }


protected:
    ObserverLb<TRandom>* M;
    Brownian<TRandom>* dW;
    TangentProcess<TRandom>* Y;
    arma::vec H;
    double a;
    bool isBs;
    double x0;
    double lambda;
    double n;
};

#endif //LOOKBACKOPTION_MARTINGALEWEIGHT_H
