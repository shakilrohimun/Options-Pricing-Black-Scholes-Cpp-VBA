#include "pch.h"
#include "OptionPricing.hpp"

// Fonction maximum en dehors des classes
double maximum(double a, double b) {
    return (a > b) ? a : b;
}

// Impl�mentation de la fonction Thomas pour r�soudre les syst�mes tridiagonaux
void thomasSolver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, std::vector<double>& d) {
    int n = d.size();
    std::vector<double> c_star(n, 0.0);
    std::vector<double> d_star(n, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    // Balayage vers l'avant
    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    // Substitution r�troactive
    for (int i = n - 2; i >= 0; i--) {
        d_star[i] -= c_star[i] * d_star[i + 1];
    }

    d = d_star;  // Solution est dans d
}

// ctor
Option::Option(double spot_, double strike_, double rate_, double volatility_, double timeToMaturity_)
    : spot(spot_), strike(strike_), rate(rate_), volatility(volatility_), timeToMaturity(timeToMaturity_) {

    std::string errorMessage;

    if (spot_ <= 0) {
        errorMessage += "Erreur : Le prix spot doit �tre un nombre positif sup�rieur � z�ro.\n";
    }
    if (strike_ <= 0) {
        errorMessage += "Erreur : Le prix d'exercice doit �tre un nombre positif sup�rieur � z�ro.\n";
    }
    if (volatility < 0) {
        errorMessage += "Erreur : La volatilit� doit �tre un nombre positif sup�rieur � z�ro.\n";
    }
    if (timeToMaturity <= 0) {
        errorMessage += "Erreur : Le temps jusqu'� l'�ch�ance doit �tre un nombre positif sup�rieur � z�ro.\n";
    }

    if (!errorMessage.empty()) {
        MessageBoxA(NULL, errorMessage.c_str(), "Validation des donn�es", MB_OK | MB_ICONERROR);
        throw std::invalid_argument("Param�tres invalides.");
    }
}

// fonction de r�partition de la loi normale standard
double Option::norm_cdf(double x) const {
    return 0.5 * std::erfc(-x * 1 / std::sqrt(2.0));
}

// densit� de probabilit� de la loi normale standard
double Option::norm_pdf(double x) const {
    constexpr double pi = 3.14159265358979323846;
    return (1.0 / std::sqrt(2.0 * pi)) * std::exp(-0.5 * x * x);
}

// formule de d1 et d2
double Option::d1() const {
    return (std::log(spot / strike) + (rate + 0.5 * std::pow(volatility, 2)) * timeToMaturity) / (volatility * std::sqrt(timeToMaturity));
}

double Option::d2() const {
    return d1() - volatility * std::sqrt(timeToMaturity);
}

double Option::get_d1() const {
    return d1();
}

double Option::get_d2() const {
    return d2();
}

// Calcul du prix de l'option
double Option::calculatePrice(const Payoff& payoff) const {
    double d1_ = d1();
    double d2_ = d2();
    if (dynamic_cast<const PayoffCall*>(&payoff)) {
        return spot * norm_cdf(d1_) - strike * std::exp(-rate * timeToMaturity) * norm_cdf(d2_);
    }
    else {
        return strike * std::exp(-rate * timeToMaturity) * norm_cdf(-d2_) - spot * norm_cdf(-d1_);
    }
}

// Calcul du Delta
double Option::calculateDelta(bool isCall) const {
    double d1_ = d1();
    return isCall ? norm_cdf(d1_) : norm_cdf(d1_) - 1.0;
}

// Calcul du Gamma
double Option::calculateGamma() const {
    double d1_ = d1();
    return norm_pdf(d1_) / (spot * volatility * std::sqrt(timeToMaturity));
}

// Calcul du Vega
double Option::calculateVega() const {
    double d1_ = d1();
    return spot * norm_pdf(d1_) * std::sqrt(timeToMaturity);
}

// Calcul du Th�ta
double Option::calculateTheta(bool isCall) const {
    double d1_ = d1();
    double d2_ = d2();
    double firstTerm = -spot * norm_pdf(d1_) * volatility / (2 * std::sqrt(timeToMaturity));
    double secondTerm = isCall ?
        -rate * strike * std::exp(-rate * timeToMaturity) * norm_cdf(d2_) :
        rate * strike * std::exp(-rate * timeToMaturity) * norm_cdf(-d2_);
    return firstTerm + secondTerm;
}

// Calcul du Rho
double Option::calculateRho(bool isCall) const {
    double d2_ = d2();
    return isCall ?
        timeToMaturity * strike * std::exp(-rate * timeToMaturity) * norm_cdf(d2_) :
        -timeToMaturity * strike * std::exp(-rate * timeToMaturity) * norm_cdf(-d2_);
}

// Nouvelle fonction de volatilit� locale
double Option::sigma(double spot) {
    return 0.2 + 0.04 * spot + 0 * spot * spot;
}

// Impl�mentation pour calculer le prix du call avec la volatilit� locale
double Option::localVolCallPrice(double spot, double strike, double rate, double maturity) {
    const int timeSteps = 100, spaceSteps = 100;
    const double dt = maturity / timeSteps, dx = 2 * strike / spaceSteps;
    std::vector<std::vector<double>> grid(timeSteps + 1, std::vector<double>(spaceSteps + 1, 0));

    // D�finir la condition de gain terminal pour une option call
    for (int j = 0; j <= spaceSteps; ++j) {
        grid[timeSteps][j] = maximum(j * dx - strike, 0.0);
    }

    // Coefficients pour le sch�ma implicite
    std::vector<double> a(spaceSteps - 1, 0.0);
    std::vector<double> b(spaceSteps - 1, 0.0);
    std::vector<double> c(spaceSteps - 1, 0.0);
    std::vector<double> d(spaceSteps - 1, 0.0);

    // M�thode des diff�rences finies (Sch�ma implicite utilisant l'algorithme de Thomas)
    for (int i = timeSteps - 1; i >= 0; --i) {
        for (int j = 1; j < spaceSteps; ++j) {
            double S = j * dx;
            double vol = sigma(spot);  // Utilisation de la nouvelle m�thode sigma(double spot)
            double alpha = 0.5 * dt * (vol * vol * S * S / (dx * dx) - rate * S / dx);
            double beta = -dt * (vol * vol * S * S / (dx * dx) + rate);
            double gamma = 0.5 * dt * (vol * vol * S * S / (dx * dx) + rate * S / dx);

            if (j == 1) {  // Conditions aux limites � la fronti�re inf�rieure
                a[j - 1] = 0.0;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = gamma;
                d[j - 1] = grid[i + 1][j] - alpha * grid[i + 1][j - 1]; // En tenant compte de la fronti�re inf�rieure
            }
            else if (j == spaceSteps - 1) {  // Conditions aux limites � la fronti�re sup�rieure
                a[j - 1] = alpha;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = 0.0;
                d[j - 1] = grid[i + 1][j] - gamma * grid[i + 1][j + 1]; // En tenant compte de la fronti�re sup�rieure
            }
            else {  // Points int�rieurs
                a[j - 1] = alpha;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = gamma;
                d[j - 1] = grid[i + 1][j];
            }
        }

        // R�soudre le syst�me tridiagonal en utilisant l'algorithme de Thomas
        thomasSolver(a, b, c, d);

        // Mettre � jour les valeurs de la grille apr�s avoir r�solu le syst�me
        for (int j = 1; j < spaceSteps; ++j) {
            grid[i][j] = d[j - 1];
        }
    }

    return std::exp(-rate * maturity) * grid[0][spaceSteps / 2];  // Prix actualis�
}

// Impl�mentation pour calculer le prix du put avec la volatilit� locale
double Option::localVolPutPrice(double spot, double strike, double rate, double maturity) {
    const int timeSteps = 100, spaceSteps = 100;
    const double dt = maturity / timeSteps, dx = 2 * strike / spaceSteps;
    std::vector<std::vector<double>> grid(timeSteps + 1, std::vector<double>(spaceSteps + 1, 0));

    // D�finir la condition de gain terminal pour une option put
    for (int j = 0; j <= spaceSteps; ++j) {
        grid[timeSteps][j] = maximum(strike - j * dx, 0.0);
    }

    // Coefficients pour le sch�ma implicite
    std::vector<double> a(spaceSteps - 1, 0.0);
    std::vector<double> b(spaceSteps - 1, 0.0);
    std::vector<double> c(spaceSteps - 1, 0.0);
    std::vector<double> d(spaceSteps - 1, 0.0);

    // M�thode des diff�rences finies (Sch�ma implicite utilisant l'algorithme de Thomas)
    for (int i = timeSteps - 1; i >= 0; --i) {
        for (int j = 1; j < spaceSteps; ++j) {
            double S = j * dx;
            double vol = sigma(spot);  // Utilisation de la nouvelle m�thode sigma(double spot)
            double alpha = 0.5 * dt * (vol * vol * S * S / (dx * dx) - rate * S / dx);
            double beta = -dt * (vol * vol * S * S / (dx * dx) + rate);
            double gamma = 0.5 * dt * (vol * vol * S * S / (dx * dx) + rate * S / dx);

            if (j == 1) {  // Conditions aux limites � la fronti�re inf�rieure
                a[j - 1] = 0.0;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = gamma;
                d[j - 1] = grid[i + 1][j] - alpha * grid[i + 1][j - 1]; // En tenant compte de la fronti�re inf�rieure
            }
            else if (j == spaceSteps - 1) {  // Conditions aux limites � la fronti�re sup�rieure
                a[j - 1] = alpha;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = 0.0;
                d[j - 1] = grid[i + 1][j] - gamma * grid[i + 1][j + 1]; // En tenant compte de la fronti�re sup�rieure
            }
            else {  // Points int�rieurs
                a[j - 1] = alpha;
                b[j - 1] = 1.0 - beta;
                c[j - 1] = gamma;
                d[j - 1] = grid[i + 1][j];
            }
        }

        // R�soudre le syst�me tridiagonal en utilisant l'algorithme de Thomas
        thomasSolver(a, b, c, d);

        // Mettre � jour les valeurs de la grille apr�s avoir r�solu le syst�me
        for (int j = 1; j < spaceSteps; ++j) {
            grid[i][j] = d[j - 1];
        }
    }

    return std::exp(-rate * maturity) * grid[0][spaceSteps / 2];  // Prix actualis�
}

// Impl�mentation des fonctions de couverture (hedging)
double Option::hedgeCall(double spot, double strike, double rate, double maturity) {
    double delta = (localVolCallPrice(spot + 0.01, strike, rate, maturity) - localVolCallPrice(spot - 0.01, strike, rate, maturity)) / 0.02;
    return delta;
}

double Option::hedgePut(double spot, double strike, double rate, double maturity) {
    double delta = (localVolPutPrice(spot + 0.01, strike, rate, maturity) - localVolPutPrice(spot - 0.01, strike, rate, maturity)) / 0.02;
    return delta;
}

// Exportation des fonctions pour VBA

// d1 et d2
extern "C" __declspec(dllexport) double get_d1(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.get_d1();
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double get_d2(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.get_d2();
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Prix Call / Put Black-Scholes
extern "C" __declspec(dllexport) double calculate_call_price(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        PayoffCall payoff(strike);
        return option.calculatePrice(payoff);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double calculate_put_price(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        PayoffPut payoff(strike);
        return option.calculatePrice(payoff);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Delta
extern "C" __declspec(dllexport) double calculate_call_delta(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateDelta(true);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double calculate_put_delta(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateDelta(false);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Gamma
extern "C" __declspec(dllexport) double calculate_gamma(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateGamma();
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Vega
extern "C" __declspec(dllexport) double calculate_vega(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateVega();
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Th�ta
extern "C" __declspec(dllexport) double calculate_call_theta(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateTheta(true);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double calculate_put_theta(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateTheta(false);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Rho
extern "C" __declspec(dllexport) double calculate_call_rho(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateRho(true);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double calculate_put_rho(double spot, double strike, double rate, double volatility, double timeToMaturity) {
    try {
        Option option(spot, strike, rate, volatility, timeToMaturity);
        return option.calculateRho(false);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Prix Call / Put Volatilit� Local
extern "C" __declspec(dllexport) double calculate_local_vol_call_price(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta) {
    try {
        Option option(spot, strike, rate, 0, timeToMaturity);
        return option.localVolCallPrice(spot, strike, rate, timeToMaturity);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double calculate_local_vol_put_price(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta) {
    try {
        Option option(spot, strike, rate, 0, timeToMaturity);
        return option.localVolPutPrice(spot, strike, rate, timeToMaturity);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

// Hedge Call / Put Volatilite local
extern "C" __declspec(dllexport) double hedge_call_option(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta) {
    try {
        Option option(spot, strike, rate, 0, timeToMaturity);
        return option.hedgeCall(spot, strike, rate, timeToMaturity);
    }
    catch (const std::exception& e) {
        return 0;
    }
}

extern "C" __declspec(dllexport) double hedge_put_option(double spot, double strike, double rate, double timeToMaturity, double sigma_0, double alpha, double beta) {
    try {
        Option option(spot, strike, rate, 0, timeToMaturity);
        return option.hedgePut(spot, strike, rate, timeToMaturity);
    }
    catch (const std::exception& e) {
        return 0;
    }
}
