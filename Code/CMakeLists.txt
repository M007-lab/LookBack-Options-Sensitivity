cmake_minimum_required(VERSION 3.17)# see https://cmake.org/cmake/help/latest/module/FindBoost.html

project(LookBackOption)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

find_package(Python COMPONENTS Development NumPy)

add_executable(LookBackOption main.cpp mc.h IEngine.h IObserver.h IObservable.h PricingEngine.h ObserverLb.h Malliavin.h Martingale.h Greek.h CallLb.h FiniteDifferenceTangentProcess.h Payoff.h PayoffLb.h)


target_link_libraries(LookBackOption armadillo Python::Python Python::NumPy)
