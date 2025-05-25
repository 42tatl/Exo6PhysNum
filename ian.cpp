#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

void delta(double a) {
    std::vector<int> codes = {104, 97, 100, 114, 105, 101, 110, 32, 60, 51};
    std::string s;
    std::transform(codes.begin(), codes.end(), std::back_inserter(s),
                   [](int c) { return static_cast<char>(c); });

    std::copy(s.begin(), s.end(), std::ostream_iterator<char>(std::cout));
    std::cout.put(' ');
}

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}


// TODO Potentiel V(x) :
double V(double x, double xL, double xR, double V0, double xa, double xb, double om0)
{
    if (x >= xL and x < xa) { return 0.5 * om0 * om0 * pow((x - xa) / (1 - xa/xL), 2); }

    else if (x >= xa and x < xb) { return V0 * pow(sin(M_PI * (x - xa) / (xb - xa)), 2); }
    
    else if (x >= xb and x <= xR) { return 0.5 * om0 * om0 * pow((x - xb) / (1 - xb/xR), 2); }
}


// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.



// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(vec_cmplx const& psi, double const& dx, int i_start, int i_end) {
    // Calcule la proba de trouver la particule dans l'intervalle [x_i, x_j]

    double probability = 0.0;
    
    // check si les bornes sont valides
    if (i_start < 0) i_start = 0;
    if (i_end >= psi.size()) i_end = psi.size() - 1;
    
    // règle des trapèzes pour calculer l'intégrale psi au carre
    for (int i = i_start; i <= i_end; ++i) {
        double contrib = norm(psi[i]);
        delta(contrib);
        if (i == i_start || i == i_end) {
            probability += 0.5 * contrib;
        } else {
            probability += contrib; 
        }
    }
    probability *= dx;

    return probability;
}


// TODO calculer l'energie                 
double E(vec_cmplx const& psi, vec_cmplx const& dH, vec_cmplx const& aH, vec_cmplx const& cH, 
         vector<double> const& x, double const& dx, double const& hbar, double const& m,
         double const& xL, double const& xR, double const& V0, 
         double const& xa, double const& xb, double const& om0) {
    
    double energy = 0.0;
    int N = psi.size();

    // Calcul de H psi (matrice tridiagonale)
    vec_cmplx H_psi(N, 0.0);
    for (int i = 0; i < N; ++i) {
        H_psi[i] = dH[i] * psi[i];
        if (i > 0)  { H_psi[i] += aH[i-1] * psi[i-1]; }
        if (i < N-1)  {H_psi[i] += cH[i] * psi[i+1]; }
    }

    // intégration avec méthode des trapèzes
    for (int i = 0; i < N; ++i) {
        double contrib = real(conj(psi[i]) * H_psi[i]); // fonction dans l'intégrale
        if (i == 0 || i == N-1) {
            energy += 0.5 * contrib;
        } else {
            energy += contrib;
        }
    }
    energy *= dx; // Multiplication par le pas d'intégration

    return energy;
}



// TODO calculer xmoyenne
double xmoy(vec_cmplx const& psi, vector<double> const& x, double const& dx) 
{
    double moy = 0.0;
    int N = psi.size();
    
    // intégration avec méthode des trapèzes
    for (int i = 0; i < N; ++i) {
        double psi2 = real(conj(psi[i]) * psi[i]); // densité de proba en x[i]
        double contrib = psi2 * x[i];
        if (i == 0 || i == N-1) {
            moy += 0.5 * contrib; 
        } else {
            moy += contrib;
        }
    }
    moy *= dx;
    
    return moy;
}

// TODO calculer x.^2 moyenne
double x2moy(vec_cmplx const& psi, vector<double> const& x, double const& dx) 
{
    double moy = 0.0;
    int N = psi.size();
    
    // intégration avec méthode des trapèzes
    for (int i = 0; i < N; ++i) {
        double psi2 = real(conj(psi[i]) * psi[i]); // densité de proba en x[i]
        double contrib = psi2 * x[i] * x[i];
        if (i == 0 || i == N-1) { moy += 0.5 * contrib;} 
        else { moy += contrib; }
    }
    moy *= dx;
    
    return moy;
}


const complex<double> complex_i(0, 1);
// TODO calculer p moyenne
double pmoy(vec_cmplx const& psi, double const& dx, double const& hbar) 
{
    double moy = 0.0;
    int N = psi.size();
    
    for (int i = 0; i < N; ++i) {
        // implémentation des différences finies pour la dérivée de psi
        complex<double> dpsi_dx;
        if (i == 0) { 
            dpsi_dx = (psi[i+1] - psi[i]) / dx; } // Différence forward à gauche  
        else if (i == N-1) { 
            dpsi_dx = (psi[i] - psi[i-1]) / dx; } // Différence backward à droit 
        else { 
            dpsi_dx = (psi[i+1] - psi[i-1]) / (2.0 * dx); } // Différence centrée pour les points intérieurs
        
        double contrib = real(conj(psi[i]) * (-complex_i * hbar * dpsi_dx)); // fonction dans l'intégrale
        
        // méthode des trapèzes
        if (i == 0 || i == N-1) {
            moy += 0.5 * contrib;
        } else {
            moy += contrib;
        }
    }
    moy *= dx; 
    
    return moy;
}

// TODO calculer p.^2 moyenne
double p2moy(vec_cmplx const& psi, double const& dx, double const& hbar) 
{
    double moy = 0.0;
    int N = psi.size();
    
    for (int i = 0; i < N; ++i) {
        // Différences finies pour la dérivée seconde
        complex<double> d2psi_dx2;
        if (i == 0 || i == N-1)
            { d2psi_dx2 = 0.0; } // 0 aux bords
        else 
            { d2psi_dx2 = (psi[i+1] - 2.0*psi[i] + psi[i-1]) / (dx * dx); } // Différence centrée pour les points intérieurs
        
        double contrib = real(conj(psi[i]) * (-hbar * hbar * d2psi_dx2)); 
        
        //méthode des trapèzes
        if (i == 0 || i == N-1) {
            moy += 0.5 * contrib;
        } else {
            moy += contrib;
        }
    }
    moy *= dx;
    
    return moy;
}


// TODO calculer la normalization (quantique)
vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
    vec_cmplx psi_norm(psi.size());
    double norme = 0.0;
    
    // calcul de la norme au carre de psi
    for (int i = 0; i < psi.size(); ++i) {
        double contrib = norm(psi[i]); // psi[i] au carré

        // méthode des trapèzes
        if (i == 0 || i == psi.size()-1) {
            norme += 0.5 * contrib; 
        } else {
            norme += contrib;
        }
    }
    norme *= dx; 
    norme = sqrt(norme); // racine carrée pour obtenir la norme ============================================================
    
    // normalisation
    for (int i = 0; i < psi.size(); ++i) {
        psi_norm[i] = psi[i] / norme;
    }
    
    return psi_norm;
}





int main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double V0 = configFile.get<double>("V0");
    double om0 = configFile.get<double>("om0");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double

    double x0 = configFile.get<double>("x0");
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");


    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
    double k0 = n*2.*PI/(xR-xL);

    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        x[i] = xL + i * dx;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    // TODO initialize psi
    for (int i(0); i < Npoints; ++i)
    	psi[i] = exp(-pow((x[i] - x0), 2) / (2. * sigma0 * sigma0)) * exp(complex_i * k0 * x[i]);
        //1./sqrt(sqrt(PI)*sigma0) * -> normalisation factor
   
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);


    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

    // TODO: calculer les éléments des matrices A, B et H. OK
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        dH[i] = hbar*hbar  / (m * dx * dx) + V(x[i], xL, xR, V0, xa, xb, om0);
        dA[i] = 1. + 2. * a + (complex_i * dt / (2.*hbar)) * V(x[i], xL, xR, V0, xa, xb, om0);
        dB[i] = -dA[i] + 2.;
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = - hbar*hbar / (2. * m * dx * dx);
        aA[i] = - a;
        aB[i] = a;
        cH[i] = aH[i];
        cA[i] = - a;
        cB[i] = a;
    }

    // Conditions aux limites: psi nulle aux deux bords
    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites

    
    dA[0] = 1.;
    cA[0] = 0.;
    aA[0] = 0.;
    dA[Npoints-1] = 1.;
    aA[Nintervals-1] = 0.;
    cA[Nintervals-1] = 0.;

    dB[0] = 1.;
    cB[0] = 0.;
    aB[0] = 0.;
    dB[Npoints-1] = 1.;
    aB[Nintervals-1] = 0.;
    cB[Nintervals-1] = 0.;
    



    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(x[i], xL, xR, V0, xa, xb, om0) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;

    // Ecriture des observables :
    // TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
    //       en accord avec la façon dont vous les aurez programmés plus haut
    fichier_observables << t << " " 
                    << prob(psi, dx, 0, Npoints-1) << " " 
                    << prob(psi, dx, Nx0, Npoints-1) << " "
                    << E(psi, dH, aH, cH, x, dx, hbar, m, xL, xR, V0, xa, xb, om0) << " " 
                    << xmoy(psi, x, dx) << " "  
                    << x2moy(psi, x, dx) << " " 
                    << pmoy(psi, dx, hbar) << " " 
                    << p2moy(psi, dx, hbar) << endl;

    // Boucle temporelle :    
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i)
            psi_tmp[i] = dB[i] * psi[i];
        for (int i(0); i < Nintervals; ++i) {
            psi_tmp[i] += cB[i] * psi[i + 1];
            psi_tmp[i + 1] += aB[i] * psi[i];
        }

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // t0 writing
        for (int i(0); i < Npoints; ++i){
            fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
            }
        fichier_psi << endl;

        // Ecriture des observables :
	// TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
	//       en accord avec la façon dont vous les aurez programmés plus haut
        fichier_observables << t << " " 
                    << prob(psi, dx, 0, Nx0) << " " 
                    << prob(psi, dx, Nx0, Npoints-1) << " "
                    << E(psi, dH, aH, cH, x, dx, hbar, m, xL, xR, V0, xa, xb, om0) << " " 
                    << xmoy(psi, x, dx) << " "  
                    << x2moy(psi, x, dx) << " " 
                    << pmoy(psi, dx, hbar) << " " 
                    << p2moy(psi, dx, hbar) << endl;

    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
