// -*- coding: utf-8 -*-
#include <random>
#include "PricingEngine.h"
#include "mc.h"
#include "CallLb.h"
#include "FiniteDifferenceTangentProcess.h"
#include "quantile.h"
#include "Antithetic.h"
#include "Malliavin.h"
#include "Martingale.h"
#include "timer.hpp"

using namespace std;

int main(int argc, char** argv) {
    timer t;
    random_device rd;
    auto seed = rd();
    mt19937_64 gen{seed};
    double r = 0.05, sigma = 0.25;
    double x0 = 100;
    double T = 1;
    unsigned N = atoi(argv[2]);
    double K = 120;
    int n = atoi(argv[1]);
    auto h = T / (double) n;
    double gymin = 0.02, gymax = 0.035;
    double dymin = 0.5, dymax = 0.8;
    std::normal_distribution<> G(0, 1);
    using TRandom = std::normal_distribution<>;
    auto X = PricingEngine<TRandom>(x0, r, sigma, h, G, n);
    PayoffLb<TRandom>* call = new CallLb<TRandom>(X, K);

    std::cout << "> Sample Size = " << N << ", # steps = " << n << "\n";

    std::cout << ">>> [Delta] WithOut Variance Reduction" << "\n";
    Greek<TRandom>* greek = new Delta_DFD<TRandom>(call,10);
    call->set_greek(greek);
    plt::figure_size(1200, 780);
    t.reset();
    std::cout << "> Decreasing Finite difference: \n\t " << monte_carlo(*call,gen,N,std::exp(-r*T),"Decreasing Finite difference") << "\n";
    std::cout << t << "\n";
    delete greek;

    X.reinitObservers();
    greek = new Delta_FD<TRandom>(call);
    call->set_greek(greek);
    t.reset();
    std::cout << "> Finite difference: \n\t " << monte_carlo(*call,gen,N,std::exp(-r*T),"Finite difference") << "\n";
    delete greek;

    X.reinitObservers();
    greek = new Delta_TP<TRandom>(call);
    call->set_greek(greek);
    t.reset();
    std::cout << "> Tangent process: \n\t "  << monte_carlo(*call,gen,N,std::exp(-r*T),"Tangent process") << "\n";
    std::cout << t << "\n\n";
    delete greek;

    X.reinitObservers();
    greek = new Delta_Martingale<TRandom>(call,K-x0,5,1);
    call->set_greek(greek);
    t.reset();
    std::cout << "> Martingale: \n\t " << monte_carlo(*call,gen,N,std::exp(-r*T),"Martingale") << "\n";
    std::cout << t << "\n\n";
    delete greek;
    
    X.reinitObservers();
    greek = new Delta_Malv<TRandom>(call,log(K/x0));
    call->set_greek(greek);
    t.reset();
    std::cout << "> Localized Malliavin : \n\t" << monte_carlo(*call,gen,N,std::exp(-r*T),"Malliavin Localized",call->delta(x0,r,sigma,T),dymin,dymax) << "\n";
    std::cout << t << "\n\n";
    delete greek;
    plt::title("Delta (without Var. Reduction)");
    plt::save("./Grph_delta_wo_vr.png");

    //////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    
    std::cout << ">>> [Gamma] WithOut Variance Reduction" << "\n";
    X.reinitObservers();
    greek = new Gamma_FD<TRandom>(call);
    call->set_greek(greek);
    plt::figure_size(1200, 780);
    t.reset();
    std::cout << "> Finite difference: \n\t " << monte_carlo(*call,gen,N,std::exp(-r*T),"Finite difference") << "\n";
    std::cout << t << "\n\n";
    delete greek;

    X.reinitObservers();
    greek = new Gamma_FD_TP<TRandom>(call);
    call->set_greek(greek);
    t.reset();
    std::cout << "> Finite difference & Tangent Process: \n\t " << monte_carlo(*call,gen,N,std::exp(-r*T),"Finite difference & Tangent Process",call->gamma(x0,r,sigma,T),gymin,gymax) << "\n";
    std::cout << t << "\n\n";
    plt::title("Gamma (without Var. Reduction)");
    plt::save("./Grph_gamma_wo_vr.png");

    //::::With Var Reduction
    std::cout << ">>> [Delta] With Variance Reduction" << "\n";
    greek = new Delta_DFD<TRandom>(call,10);
    call->set_greek(greek);
    auto anti_call = new Antithetic<TRandom>(call);
    plt::figure_size(1200, 780);
    t.reset();
    std::cout << "> Decreasing Finite difference: \n\t" << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Decreasing Finite difference") << "\n";
    std::cout << t << "\n\n";
    delete anti_call;
    delete greek;

    X.reinitObservers();
    greek = new Delta_FD<TRandom>(call);
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    t.reset();
    std::cout << "> Finite difference: \n\t" << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Finite difference") << "\n";
    delete anti_call;
    delete greek;

    X.reinitObservers();
    greek = new Delta_TP<TRandom>(call);
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    t.reset();
    std::cout << "> Tangent process: \n\t "  << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Tangent process") << "\n";
    std::cout << t << "\n\n";
    delete anti_call;
    delete greek;

    X.reinitObservers();
    greek = new Delta_Martingale<TRandom>(call,K-x0,5,1);
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    std::cout << "> Martingale: \n\t " << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Martingale") << "\n";
    delete anti_call;
    delete greek;

    X.reinitObservers();
    greek = new Delta_Malv<TRandom>(call,log(K/x0));
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    t.reset();
    std::cout << "> Localized Malliavin: \n\t " << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Malliavin Localized",call->delta(x0,r,sigma,T),dymin,dymax) << "\n";
    std::cout << t << "\n\n";
    delete anti_call;
    delete greek;
    plt::title("Delta (with Var. Reduction)");
    plt::save("./Grph_delta_with_vr.png");


    std::cout << ">>> [Gamma] With Variance Reduction" << "\n";
    X.reinitObservers();
    greek = new Gamma_FD<TRandom>(call);
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    plt::figure_size(1200, 780);
    t.reset();
    std::cout << "> Finite difference: \n\t" << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Finite difference") << "\n";
    std::cout << t << "\n\n";
    delete anti_call;
    delete greek;

    X.reinitObservers();
    greek = new Gamma_FD_TP<TRandom>(call);
    call->set_greek(greek);
    anti_call = new Antithetic<TRandom>(call);
    t.reset();
    std::cout << "> Finite difference & Tangent Process: \n\t " << monte_carlo(*anti_call,gen,N,std::exp(-r*T),"Finite difference & Tangent Process",call->gamma(x0,r,sigma,T),gymin,gymax) << "\n";
    std::cout << t << "\n\n";
    plt::title("Gamma (with Var. Reduction)");
    plt::save("./Grph_gamma_with_vr.png");

    std::cout << "- Price : " << call->price(x0,r,sigma,T) << "\n";
    std::cout << "- Delta : " << call->delta(x0,r,sigma,T)<< "\n";
    std::cout << "- Gamma : " << call->gamma(x0,r,sigma,T)<< "\n";

    return 0;
};