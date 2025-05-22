#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double> > vec_cmplx;

double PI=3.1415926535897932384626433832795028841971e0;

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
double V(double xR, double xL, double xa, double xb, double V0, double om0, double x)
{
    const double eps = 1e-12;

    // Cas harmonique pur sur tout le domaine
    if (V0 == 0 && fabs(xa) < eps && fabs(xb) < eps) {
        if (x < 0) {
            return 0.5 * om0 * om0 * pow(x / (xL), 2);
        } else {
            return 0.5 * om0 * om0 * pow(x / (xR), 2);
        }
    }

    
    if (xL <= x && xa > x) {
        return 0.5*om0*om0*((x-xa)/(xL-xa))*((x-xa)/(xL-xa));
    } else if (x >= xa && x <= xb) {
        return V0 * sin(PI*(x - xa) / (xb - xa))*sin(PI*(x - xa) / (xb - xa));
    } else if (xb < x && xR >= x) {
        return 0.5*om0*om0*((x-xb)/(xR-xb))*((x-xb)/(xR-xb));
    }
    else
        return 10e10; // Potentiel infini

    
}


// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(double x_i, double x_j, vec_cmplx psi,vector<double> const& x, double const& dx)
{
    double probability = 0.0;
    // Parcourir les points de maillage
    for (size_t i = 0; i < x.size() - 1; ++i) {
        // Vérification  si le segment [x[i], x[i+1]] est dans l'intervalle [x_i, x_j]
        if (x[i] >= x_i && x[i + 1] <= x_j) {
            // Méthode des trapèzes : (|psi[i]|^2 + |psi[i+1]|^2) * dx / 2
            probability += (std::norm(psi[i]) + std::norm(psi[i + 1])) * dx / 2.0;
        }
    }

    return probability;
}

// TODO calculer l'energie
double E(vec_cmplx const& psi, vec_cmplx const& dH, vec_cmplx const& aH, vec_cmplx const& cH, double dx)
{
    size_t N = psi.size();
    vec_cmplx H_psi(N, 0.0); // Vecteur pour stocker H * psi

    // Appliquer H à psi
    for (size_t i = 0; i < N; ++i) {
        H_psi[i] = dH[i] * psi[i]; // Contribution de la diagonale principale
        if (i > 0) {
            H_psi[i] += aH[i - 1] * psi[i - 1]; // Contribution de la sous-diagonale
        }
        if (i < N - 1) {
            H_psi[i] += cH[i] * psi[i + 1]; // Contribution de la sur-diagonale
        }
    }

    // Calcul de l'énergie moyenne avec la méthode des trapèzes
    double energy = 0.0;
    for (size_t i = 0; i < N - 1; ++i) {
        // Méthode des trapèzes : (psi[i]* conj(H_psi[i]) + psi[i+1]* conj(H_psi[i+1])) * dx / 2
        energy += 0.5 * dx * (real(conj(psi[i]) * H_psi[i]) + real(conj(psi[i + 1]) * H_psi[i + 1]));
    }
    //cout<<"energy = " << energy << endl;
    return energy;
}

// TODO calculer xmoyenne
double xmoy(vec_cmplx const& psi, vector<double> const& x, double dx)
{
    size_t N = psi.size();
    double xmoyenne = 0.0;

    // Calcul de la position moyenne avec la méthode des trapèzes
    for (size_t i = 0; i < N - 1; ++i) {
        // Méthode des trapèzes : (x[i] * |psi[i]|^2 + x[i+1] * |psi[i+1]|^2) * dx / 2
        xmoyenne += 0.5 * dx * (x[i] * std::norm(psi[i]) + x[i + 1] * std::norm(psi[i + 1]));
    }

    return xmoyenne;
}

// TODO calculer x.^2 moyenne
double x2moy(vec_cmplx const& psi, vector<double> const& x, double dx)
{
    size_t N = psi.size();
    double x2moyenne = 0.0;

    // Calcul de la position moyenne avec la méthode des trapèzes
    for (size_t i = 0; i < N - 1; ++i) {
        // Méthode des trapèzes : (x[i] * |psi[i]|^2 + x[i+1] * |psi[i+1]|^2) * dx / 2
        x2moyenne += 0.5 * dx * (x[i]*x[i] * std::norm(psi[i]) + x[i + 1]*x[i+1] * std::norm(psi[i + 1]));
    }

    return x2moyenne;
}

// TODO calculer p moyenne
double pmoy(vec_cmplx const& psi, double dx, double hbar)
{
    size_t N = psi.size();
    vec_cmplx dpsi_dx(N, 0.0); // Vecteur pour stocker la dérivée spatiale de psi

    // Calcul de la dérivée spatiale de psi
    for (size_t i = 1; i < N - 1; ++i) {
        dpsi_dx[i] = (psi[i + 1] - psi[i - 1]) / (2.0 * dx); // Différences centrées
    }
    dpsi_dx[0] = (psi[1] - psi[0]) / dx; // Différences forward pour le bord gauche
    dpsi_dx[N - 1] = (psi[N - 1] - psi[N - 2]) / dx; // Différences backward pour le bord droit

    // Calcul de la quantité de mouvement moyenne avec la méthode des trapèzes
    double pmoyenne = 0.0;
    for (size_t i = 0; i < N - 1; ++i) {
        // Contribution de chaque segment avec la méthode des trapèzes
        pmoyenne += 0.5 * dx * (
            real(conj(psi[i]) * (-complex<double>(0, 1) * hbar * dpsi_dx[i])) +
            real(conj(psi[i + 1]) * (-complex<double>(0, 1) * hbar * dpsi_dx[i + 1]))
        );
    }

    return pmoyenne;
}

// TODO calculer p.^2 moyenne
double p2moy(vec_cmplx const& psi, double dx, double hbar)
{
    size_t N = psi.size();
    vec_cmplx d2psi_dx2(N, 0.0); // Vecteur pour stocker la dérivée seconde de psi

    // Calcul de la dérivée seconde de psi
    for (size_t i = 1; i < N - 1; ++i) {
        d2psi_dx2[i] = (psi[i + 1] - 2.0 * psi[i] + psi[i - 1]) / (dx * dx); // Différences centrées
    }
    d2psi_dx2[0] = 0.0; // Bord gauche : dérivée seconde nulle
    d2psi_dx2[N - 1] = 0.0; // Bord droit : dérivée seconde nulle

    // Calcul de p^2 moyen avec la méthode des trapèzes
    double p2moyenne = 0.0;
    for (size_t i = 0; i < N - 1; ++i) {
        // Contribution de chaque segment avec la méthode des trapèzes
        p2moyenne += 0.5 * dx * (
            real(conj(psi[i]) * (-hbar * hbar * d2psi_dx2[i])) +
            real(conj(psi[i + 1]) * (-hbar * hbar * d2psi_dx2[i + 1]))
        );
    }

    return p2moyenne;
}

// TODO calculer la normalization
vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{   


    size_t N = psi.size();
    double norm = 0.0;

    // Intégration avec la règle des trapèzes
    norm += std::norm(psi[0]);
    for (size_t i = 1; i < N - 1; ++i)
        norm += 2.0 * std::norm(psi[i]);
    norm += std::norm(psi[N - 1]);

    norm *= dx / 2.0;
    cout<<"norm = " << norm << endl;
    // computing the constant of normailzation 
    double C = 1.0 / sqrt(norm);

    // initialization of the output vector 
    vec_cmplx psi_norm(N);
    //computing the output vector
    for (size_t i = 0; i < N; ++i)
        psi_norm[i] = C * psi[i];
    

    return psi_norm;
}

//function to compute l'incertitude de p et x 
double delta_x(vec_cmplx const& psi, double dx, double hbar, vector<double> const& x)
{
    double x2moyenne= x2moy(psi, x, dx);
    double xmoyenne= xmoy(psi, x, dx);
    return sqrt(x2moyenne - xmoyenne * xmoyenne);
}

double delta_p(vec_cmplx const& psi, double dx, double hbar)
{
    double p2moyenne = p2moy(psi, dx, hbar);
    double pmoyenne = pmoy(psi, dx, hbar);
    return sqrt(p2moyenne - pmoyenne * pmoyenne);
}



int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in"); // Fichier d'input par defaut
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
    double k0 = n * PI / (xR - xL); // k0 = n * pi / L

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

    	psi[i] = exp(complex<double>(0.0, k0 * x[i]))*exp(-(x[i]-x0)*(x[i]-x0)/(2.*sigma0*sigma0));
   
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);
    for (size_t i = 0; i < psi.size(); ++i) {
        if (std::isnan(real(psi[i])) || std::isnan(imag(psi[i]))) {
            std::cerr << "💥 psi[" << i << "] est NaN juste après initialisation." << std::endl;
            exit(1);
        }
    }

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {   
        complex <double> b = complex_i  * dt * V(xR,xL,xa,xb,V0,om0,x[i])/ (2.* hbar);
       
        dH[i] = 1.0 / (m * dx * dx) + V(xR,xL,xa,xb,V0,om0,x[i]);
        dA[i] = 1.0 +2.0*a + b;
        dB[i] = 1.0 - 2.0*a - b;
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = - 1.0/(2.0*m*dx*dx); // hbar = 1;
        aA[i] = -a;
        aB[i] = a;
        cH[i] = - 1.0/(2.0*m*dx*dx); // hbar = 1
        cA[i] = -a;
        cB[i] = a;
    }

    // Conditions aux limites: psi nulle aux deux bords
    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
    dA[0] = 1.;
    cA[0] = 0.;
    aA[0] = 0.;
    dA[Npoints - 1] = 1.;
    aA[Npoints - 2] = 0.;
    cA[Npoints - 2] = 0.;

    dB[0] = 1.;
    cB[0] = 0.;
    aB[0] = 0.;
    dB[Npoints - 1] = 1.;
    aB[Npoints - 2] = 0.;
    cB[Npoints - 2] = 0.;






    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " <<V(xR,xL,xa,xb,V0,om0,x[i]) << endl;
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
    fichier_observables << t << " " << prob(xL, 0, psi, x, dx) << " " << prob(0, xR, psi, x, dx)
                << " " << E(psi,  dH,  aH,  cH, dx) << " " << xmoy (psi, x, dx) << " "  
                << x2moy(psi, x, dx) << " " << pmoy (psi, dx, hbar) << " " << p2moy(psi, dx, hbar)<< " " << delta_x( psi,  dx, hbar, x) << " " << delta_p( psi,  dx, hbar)<<endl; 

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
        fichier_observables << t << " " << prob(xL, 0, psi, x, dx) << " " << prob(0, xR, psi, x, dx)
                    << " " << E(psi,  dH,  aH,  cH, dx) << " " << xmoy (psi, x, dx) << " "  
                    << x2moy(psi, x,dx) << " " << pmoy (psi, dx, hbar) << " " << p2moy(psi, dx, hbar)<< " " << delta_x( psi,  dx, hbar, x) << " " << delta_p( psi,  dx, hbar) << endl; 
        
    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}