#ifndef OPTION_PRICING_HPP
#define OPTION_PRICING_HPP

#include <cmath>        // utile pour l'utilisation de fonction (log, sqrt, exp, erfc)
#include <vector>       // utile pour localVolPrice et thomasSolver
#include <stdexcept>    // utile pour les exportation d'erreurs
#include <string>       // utile pour générer qu'un message d'erreur

// Déclaration de la fonction maximum en dehors de toute classe pour pouvoir l'utiliser dans les classes PayoffCall/Put
double maximum(double a, double b);

// Déclaration de la fonction Thomas Solver
void thomasSolver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, std::vector<double>& d);

class Payoff {
public:
    virtual ~Payoff() {}
    virtual double operator()(double spot) const = 0;
};

class PayoffCall : public Payoff {
private:
    double strike;
public:
    PayoffCall(double strike_) : strike(strike_) {}
    virtual double operator()(double spot) const override {
        return maximum(spot - strike, 0.0);
    }
};

class PayoffPut : public Payoff {
private:
    double strike;
public:
    PayoffPut(double strike_) : strike(strike_) {}
    virtual double operator()(double spot) const override {
        return maximum(strike - spot, 0.0);
    }
};

class Option {
private:
    double spot;
    double strike;
    double rate;
    double volatility;
    double timeToMaturity;

    double norm_cdf(double x) const;
    double norm_pdf(double x) const;
    double d1() const;
    double d2() const;

public:
    Option(double spot_, double strike_, double rate_, double volatility_, double timeToMaturity_);

    double calculatePrice(const Payoff& payoff) const;

    double get_d1() const;
    double get_d2() const;

    double calculateDelta(bool isCall) const;
    double calculateGamma() const;
    double calculateVega() const;
    double calculateTheta(bool isCall) const;
    double calculateRho(bool isCall) const;

    double sigma(double spot);

    double localVolCallPrice(double spot, double strike, double rate, double maturity);
    double localVolPutPrice(double spot, double strike, double rate, double maturity);

    double hedgeCall(double spot, double strike, double rate, double maturity);
    double hedgePut(double spot, double strike, double rate, double maturity);
};

// Prototypes pour les fonctions exportées

// d1 et d2
extern "C" __declspec(dllexport) double get_d1(double spot, double strike, double rate, double volatility, double timeToMaturity);
extern "C" __declspec(dllexport) double get_d2(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Prix Call / Put Black-Scholes
extern "C" __declspec(dllexport) double calculate_call_price(double spot, double strike, double rate, double volatility, double timeToMaturity);
extern "C" __declspec(dllexport) double calculate_put_price(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Delta
extern "C" __declspec(dllexport) double calculate_call_delta(double spot, double strike, double rate, double volatility, double timeToMaturity);
extern "C" __declspec(dllexport) double calculate_put_delta(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Gamma
extern "C" __declspec(dllexport) double calculate_gamma(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Vega
extern "C" __declspec(dllexport) double calculate_vega(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Theta
extern "C" __declspec(dllexport) double calculate_call_theta(double spot, double strike, double rate, double volatility, double timeToMaturity);
extern "C" __declspec(dllexport) double calculate_put_theta(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Rho
extern "C" __declspec(dllexport) double calculate_call_rho(double spot, double strike, double rate, double volatility, double timeToMaturity);
extern "C" __declspec(dllexport) double calculate_put_rho(double spot, double strike, double rate, double volatility, double timeToMaturity);

// Prix Call / Put Volatilité local
extern "C" __declspec(dllexport) double calculate_local_vol_call_price(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta);
extern "C" __declspec(dllexport) double calculate_local_vol_put_price(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta);

// Hedge Volatilite local
extern "C" __declspec(dllexport) double hedge_call_option(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta);
extern "C" __declspec(dllexport) double hedge_put_option(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta);

#endif
