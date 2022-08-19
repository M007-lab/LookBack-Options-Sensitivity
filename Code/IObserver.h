#ifndef LOOKBACKOPTION_IOBSERVER_H
#define LOOKBACKOPTION_IOBSERVER_H

#include <armadillo>
#include <cmath>

class IObserver {
public:
    IObserver(){};
    IObserver(int n);
    IObserver(const IObserver &Observer);
    virtual void reinit();
    virtual void update(double z,double u) = 0;
    virtual double operator[](int id) const;
    virtual arma::vec get_path() const;
protected:
    arma::vec path;
    int i = 0; // current_state
};


IObserver::IObserver(int n){
    path = arma::vec(n+1);
}

IObserver::IObserver(const IObserver &Observer){
    path = Observer.get_path();
    path *= 0;
}

arma::vec IObserver::get_path() const {
    return path;
}

void IObserver::reinit(){
    double x0 = path[0];
    path = 0*path;
    i = 0;
    path[0] = x0;
}

double IObserver::operator[](int id) const {return this->path[id];}
#endif //LOOKBACKOPTION_IOBSERVER_H
