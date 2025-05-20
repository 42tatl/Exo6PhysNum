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

// Potentiel V(x) : @TODO write potential
double V(double x, double V0, double n_v, double xL, double xR, string pot)
{
    if (pot=="harm"){
        return V0*pow(2*x/(xR-xL),2);
    }
    else if (pot=="wall"){
        if (x<0) return 0;
        else return V0;
    }else if (pot=="box"){
        if (x>= xL and x<=xR) return 0.0;
        else return 1000000000;
    }
    else return 0.5*V0*(1.0+cos(2.0*M_PI*n_v*(x-xL)/(xR-xL)));
}

// @TODO compute the folliwing quantities
// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule entre les points nL.dx et nR.dx,
//  - E:    calcule son energie moyenne,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx,double);

// Les definitions de ces fonctions sont en dessous du main.

double prob(vec_cmplx psi, double bord, int Nx0,double dx,vector<double> x)
{
// pas du tout sur
    double P(0.);
    if (bord<x[Nx0]) {
        for (int i = 0; i < Nx0; ++i) {
            P += dx * 0.5 * (norm(psi[i])  + norm(psi[i + 1]) );
        }
    }else {
        for (int i = Nx0; i < psi.size()-1; ++i) {
            P += dx * 0.5 * (norm(psi[i])  + norm(psi[i + 1]));
        }
    }

    return P;
}

double E(vec_cmplx psi, vec_cmplx dH, vec_cmplx aH,vec_cmplx cH,  double dx)
{
    vec_cmplx Hpsi = vec_cmplx(psi.size(),0.0);
    for (int i= 0; i<psi.size();++i){
        Hpsi[i]= dH[i]*psi[i];
        if (i>0){
            Hpsi[i] += aH[i-1]*psi[i-1];
        } if (i<psi.size()-1) {
            Hpsi[i] += cH[i+1]*psi[i+1];
        }

    }
    double E(0.);
    for (int i= 0; i<psi.size()-1;++i){
        E+= 0.5*dx*real((conj(psi[i])*Hpsi[i] + conj(psi[i+1])*Hpsi[i+1]));
    }
    return E;
}

double xmoy(vec_cmplx psi, vector<double> x, double dx)
{
    double xm(0.);
    for (int i= 0; i<psi.size();++i){
        xm+= real(dx*0.5*(conj(psi[i])*x[i]*psi[i] +conj(psi[i+1])*x[i+1]*psi[i+1]));
    }
    return xm;
}

double x2moy(vec_cmplx psi, vector<double> x, double dx)
{
    double xm2(0.);
    for (int i= 0; i<psi.size();++i){
        xm2+= real(dx*0.5*(conj(psi[i])*pow(x[i],2)*psi[i] +conj(psi[i+1])*pow(x[i+1],2)*psi[i+1]));
    }
    return xm2;
}

double pmoy(vec_cmplx psi, double const& hbar)
{
    double pm(0.);
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    for (int i= 1; i<psi.size()-2;++i){
        pm+= real(0.25*(conj(psi[i])*(-complex_i*hbar) *(psi[i+1]-psi[i-1]) +conj(psi[i+1])*(-complex_i*hbar) *(psi[i+2]-psi[i])));
    }
    int l=psi.size();
    pm+=real(0.25*(conj(psi[0])*(-complex_i*hbar) *(psi[1]-psi[0]) +conj(psi[1])*(-complex_i*hbar) *(psi[2]-psi[0])));
    pm+=real(0.25*(-complex_i*hbar)*(conj(psi[l-1])*(psi[l] - psi[l-2] ) + 2.0*conj(psi[l])*(psi[l] - psi[l-1]) ));
    return pm;
}

double p2moy(const vec_cmplx& psi, double dx, double hbar)
{
    size_t N = psi.size();
    double result = 0.0;
    double dx2 = dx * dx;

    for (size_t i = 1; i < N - 1; ++i) {
        complex<double> d2psi = (psi[i + 1] - 2.0 * psi[i] + psi[i - 1]) / dx2;
        result += real(conj(psi[i]) * (-hbar * hbar * d2psi));
    }

    return dx * result;
}

//@TODO write a function to normalize psi

vec_cmplx normalize(vec_cmplx psi, double dx)
{
    // TODO: Normalisation de la fonction d'onde initiale psi

    vec_cmplx psi_norm(psi.size(), 0.);
    double C(0.0);
    for(int i(0);i<psi.size()-1;++i){
        //C+=0.5*dx*real(conj(psi[i])*psi[i] + conj(psi[i+1])*psi[i+1]  );
        C+=0.5*dx*(norm(psi[i]) + norm(psi[i+1]) );
    }
    for(int i(0);i<psi.size();++i){
        psi_norm[i]=psi[i]/sqrt(C);
    }
    return psi_norm;
}
int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    complex<double> complex_1= complex<double>(1.0,0.0);

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
    double V0 = configFile.get<double>("V0");
    double n_v = configFile.get<double>("n_v");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double
    string pot = configFile.get<string>("pot");
    double x_da = configFile.get<double>("xda");
    double x_db  = configFile.get<double>("xdb"); //
    string detec = configFile.get<string>("detec");
    double t_detec  = configFile.get<double>("t_detec"); // Read mode number as integer, convert to double
    string init = configFile.get<string>("init");
    // Parametres numeriques :

    double dt = configFile.get<double>("dt");
    int Nintervals = configFile.get<int>("Nintervals");
    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    //@TODO build the x mesh
    for (int i = 0; i < Npoints; ++i){
        x[i] = xL + i*dx;
    }
    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0=round((Npoints-1)/2); //chosen xR*0.5 since top of potential is at half x domain
  
    double x0 = configFile.get<double>("x0");
    double k0 = 2 * M_PI * n / (xR - xL);
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);
    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
    if (init=="Gauss"){
        for (int j = 0; j < Npoints; ++j){
            psi[j] = exp(complex_i*k0*x[j])*exp(-0.5*pow((x[j]-x0)/sigma0,2));
        }
    } else if (init=="box"){
        for (int j = 0; j < Npoints; ++j){
            int nn=n;
            if (nn%2==0) psi[j] = sin(n*M_PI*x[j]/(xR-xL)); //fct propre pour particule dans puit
            else psi[j]=cos(n*M_PI*x[j]/(xR-xL));

        }
    }

    // TODO: Modifications des valeurs aux bords :
    psi.front()=0.0;
    psi.back()=0.0;
    // TODO Normalisation :
    psi = normalize(psi,dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales supérieures et inférieures

    for (int i = 0; i < Npoints; ++i){  // pour les diagonales
        dH[i]=hbar*hbar/(m*dx*dx) + V(x[i],V0,n_v,xL,xR,pot);
        dA[i]= complex_1+ 2.0*a + complex_i*dt*V(x[i],V0,n_v,xL,xR,pot)/(hbar*2.);
        dB[i]=complex_1 - 2.0*a - complex_i*dt*V(x[i],V0,n_v,xL,xR,pot)/(hbar*2.);
    }
    for (int i = 0; i < Nintervals; ++i){ //pour les sur/sous diagonales
        aH[i]=-hbar*hbar/(2.0*m*dx*dx);
        cH[i]=-hbar*hbar/(2.0*m*dx*dx);
        aA[i]=-a;
        cA[i]=-a;
        aB[i]=a;
        cB[i]=a;
    }

    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
    aA.front()=0.0;
    dA.front()=complex_1;
    cA.front()=0.0;
    aA.back()=0.0;
    dA.back()=complex_1;
    cA.back()=0.0;

    aB.front()=0.0;
    dB.front()=complex_1;
    cB.front()=0.0;
    aB.back()=0.0;
    dB.back()=complex_1;
    cB.back()=0.0;


    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(x[i],V0,n_v,xL,xR,pot) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi2((output + "_psi2").c_str());
    fichier_psi2.precision(15);
    ofstream fichier_Repsi((output + "_Repsi").c_str());
    fichier_Repsi.precision(15);
    ofstream fichier_Impsi((output + "_Impsi").c_str());
    fichier_Impsi.precision(15);
    ofstream fichier_observables((output + "_obs").c_str());
    fichier_observables.precision(15);


    // t0 writing
    /*
    for (int i(0); i < Npoints; ++i){
        fichier_psi << pow(abs(psi[i]), 2)  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;
     */

    for (int i(0); i < Npoints; ++i){
        fichier_psi2 << norm(psi[i])  << " ";
    }
    fichier_psi2 << endl;
    for (int i(0); i < Npoints; ++i){
        fichier_Repsi << real(psi[i])  << " ";
    }
    fichier_Repsi << endl;
    for (int i(0); i < Npoints; ++i){
        fichier_Impsi << imag(psi[i])  << " ";
    }
    fichier_Impsi << endl;

    double xc=(xR+xL)/2.0;
    // Ecriture des observables :
    fichier_observables << t << " " << prob(psi,xL,Npoints,dx,x) << " " << prob(psi,xR,Npoints,dx,x) << " " << E(psi,dH,aH,cH,dx) << " " << xmoy(psi,x,dx) << " "
                << x2moy(psi,x,dx) << " " << pmoy (psi,hbar) << " " << p2moy(psi,dx,hbar) << endl;

    // Boucle temporelle :

    while (t < tfin) {

        // TODO Calcul du membre de droite :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i= 0; i<psi.size();++i){
            psi_tmp[i]+=dB[i]*psi[i];//+aB[i-1]*psi[i-1]+cB[i+1]*psi[i+1];
            if(i>0) psi_tmp[i]+=aB[i-1]*psi[i-1];
            if(i<psi_tmp.size()-1) psi_tmp[i]+=cB[i]*psi[i+1];
        }
        //cout<<t<<endl;
        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        //dt=min(dt,tfin-0.5*dt-t);
        t += dt;
        if(detec=="active"){
            if (t<t_detec+dt/2 and t >t_detec-dt/2 ){ //(t-dt/2>t_detec and t+dt <t_detec) and prob(psi,xR,Nx0,dx,x)>0.0
                cout<<"entre dans la boucle detecteur"<<endl;
                for (int i=0;i<x.size();++i){
                    if (x[i]<0){
                        psi[i]=0.0;
                    } if (x[i] >= 0.0 and x[i] < x_da){
                        psi[i] *= pow(sin(M_PI*x[i]/(2*x_da)),2);
                    } if (x[i] >= x_da and x[i] < x_db){
                        psi[i]*= 1.0;
                    } if (x[i] >= x_db and x[i] < xR){
                        psi[i]*= pow(cos(M_PI*(x[i]-x_db)/(2*(xR-x_db))),2);
                    }
                }
                psi = normalize(psi,dx);
            }
        }

        // t0 writing

        for (int i(0); i < Npoints; ++i){
            fichier_psi2 << norm(psi[i])  << " ";
        }
        fichier_psi2 << endl;
        for (int i(0); i < Npoints; ++i){
            fichier_Repsi << real(psi[i])  << " ";
        }
        fichier_Repsi << endl;
        for (int i(0); i < Npoints; ++i){
            fichier_Impsi << imag(psi[i])  << " ";
        }
        fichier_Impsi << endl;

        // Ecriture des observables :
        fichier_observables << t << " " << prob(psi,xL,Nx0,dx,x) << " " << prob(psi,xR,Nx0,dx,x) << " " << E(psi,dH,aH,cH,dx) << " " << xmoy(psi,x,dx) << " "
                            << x2moy(psi,x,dx) << " " << pmoy (psi,hbar) << " " << p2moy(psi,dx,hbar) << endl;

    } // Fin de la boucle temporelle



    fichier_observables.close();
    fichier_Repsi.close();
    fichier_Impsi.close();
    fichier_psi2.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
