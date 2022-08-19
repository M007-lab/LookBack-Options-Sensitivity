#include "list"
#include "IObserver.h"

struct IObservable{
    virtual void attach(IObserver *observer);
    virtual void detach(IObserver *observer);
    virtual void notify(double z,double u);
protected :
    std::list<IObserver *> observers;
};

void IObservable::attach(IObserver *observer){this->observers.push_back(observer);}

void IObservable::detach(IObserver *observer) {this->observers.remove(observer);}

void IObservable::notify(double z, double u){
    auto iterator = observers.begin();
    while (iterator != observers.end()) {
        (*iterator)->update(z,u);
        ++iterator;
    }
}