//
//  blitztest.cpp
//
//
//  Created by Jorge Antonio Perez Hernandez on 5/13/14.
//
//

#define __restrict__

//#define PARALLEL
#define FILE0
#define OBLATE

#define numberOfPeriods 1000
#define numberOfInitialConditions 1000

#ifdef PARALLEL
#include "mpi.h"
#endif

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <blitz/array.h>

using namespace blitz;

typedef Array<double, 3> blitz3Djet;
typedef Array<double, 2> blitz2Djet;
typedef Array<double, 1> blitz1Djet;

//Data taken from http://ssd.jpl.nasa.gov
//Saturn GM value (km^3s^-2) is M_Sa = 37,931,207.7
#define M_Su  3498.754931 //(http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
#define M_Sa  1.0
#define M_Pr  2.83144e-10 //Prometheus mass (in Saturn masses) calculated from data in http://ssd.jpl.nasa.gov/?sat_phys_par
#define M_Pa  2.43599e-10 //Pandora  mass (in Saturn masses) calculated from data in http://ssd.jpl.nasa.gov/?sat_phys_par
#define M_Ti  2.38827e-04 //Titan  mass (in Saturn masses) calculated from data in http://ssd.jpl.nasa.gov/?sat_phys_par

#define a_Sa 23785.922878 //Saturn semimajor axis wrt Sun (http://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html)
#define a_Ti 20.27386 // Titan semimajor axis (ssd.jpl.nasa.gov)
#define a_Pr 2.313483 // Prometheus's semimajor axis
#define a_Pa 2.352990 // Pandora's semimajor axis
#define a_F 2.320070 // F ring's axis
#define e_Sa 0.0565 // Saturn eccentricity wrt Sun (http://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html)
#define e_Ti 0.0288           // Titan eccentricity
#define e_Pr 0.0022          // Prometheus eccentricity
#define e_Pa 0.0042          // Pandora eccentricity
#define e_F 0.0026          // F Ring eccentricity

#define R_Hill_Pr 0.001053 //=27.712714/R_Ur (Km)
#define R_Hill_Pa 0.001019 //=31.827807/R_Ur (Km)

#define a_RP_min a_Pr //a_Pr-0.01;
#define a_RP_MAX a_Pa //a_Pa+0.01;
#define e_RP_min 0.0
#define e_RP_MAX 0.01

//number of MASSIVE particles, including test particle it's +1
#define numberOfParticles 5

#define absoluteErrorDefinition 1.0e-20
#define relativeErrorDefinition 1.0e-20
#define safetyFactor -7.5
#define timeStepControlMethodDefinition 2
#define timeStep 0.5
#define maxOrder 40
#define regulaFalsiTolerance 1.0e-14
#define newtonRaphsonTolerance 1.0e-14

#define myIndex 5

#define index_su 0
#define index_sa 1
#define index_ti 2
#define index_pr 3
#define index_pa 4
#define index_rp 5

#define rho_i index_sa
#define rho_j index_rp

#define pi   3.141592653589793  //=atan(1.0)*4.0 Let pi=1 :P
#define T_Pr 6.283185307179586  //Prometheus's period =2*pi

#ifdef OBLATE
    #define J2   0.01629071        //Jacobson, 2006, Gravity field of Saturnian System...
    #define Lambda 0.024436065     //=1.5*J2
    #define G   12.3259565 // <-with J2 //=a_Pr*a_Pr*a_Pr/((1.+M_Pr)*(1.+Lambda/(a_Pr*a_Pr))) //Universal Gravitation Constant
#else
	#define J2   0.0              //spherical planet!
	#define Lambda 0.0            //=1.5*J2
    #define G    12.382234 // <-without J2 //=a_Pr*a_Pr*a_Pr/(1.+M_Pr) //Universal Gravitation Constant
#endif

#include "PRNplus1BP_blitz.h"

int main () {

    std::cout << std::setprecision(16);

#ifndef PARALLEL
	int rank = 0; int numtasks = 1;
#endif

#ifdef PARALLEL
	int rank;
    int  numtasks, namelen, rc;
	rc = MPI_Init(NULL,NULL);
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	//Checking succesful MPI initialization...
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	//Starting MPI...
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Get_processor_name(processor_name,&namelen);
#endif

    int maxperiods=numberOfPeriods;

    int snapshots1=ceil(1.*maxperiods/100.);
    int snapshots2=ceil(5.*maxperiods/100.);
    int snapshots3=ceil(1.*maxperiods/10.);
    int snapshots4=ceil(5.*maxperiods/10.);

    if(rank==0){
        //std::cout << "rank = 0 :" << std::endl;
		std::cout << "This is a Taylor method, MPI-parallelized, PR(N+1)BP simulation." << std::endl;
        std::cout << "Integrating " << numberOfInitialConditions << " orbits; T_Max/T_Pr~" << maxperiods << "." << std::endl;
		std::cout << "Number of tasks: "<< numtasks <<"."<<std::endl;
		std::cout << "macheps=std::numeric_limits<double>::epsilon()=" << std::numeric_limits<double>::epsilon() << std::endl;
		std::cout << "J2    =" << J2     << std::endl;
		std::cout << "Lambda=" << Lambda << std::endl;
		std::cout << "G     =" << G      << std::endl;
        std::cout << "Gval  =" << a_Pr*a_Pr*a_Pr/((1.+M_Pr)*(1.+Lambda/(a_Pr*a_Pr))) << std::endl;
        //std::cout << "T_Pr=" << T_Pr << std::endl;
        std::cout << "snap1=" << snapshots1 << std::endl;
        std::cout << "snap2=" << snapshots2 << std::endl;
        std::cout << "snap3=" << snapshots3 << std::endl;
        std::cout << "snap4=" << snapshots4 << std::endl;
        std::cout << "snap5=" << numberOfPeriods << std::endl;

	};

#ifdef PARALLEL
    std::cout << "This is process #" << rank << " on " << processor_name << "." << std::endl;
#endif

    int P = numberOfInitialConditions;
    int N = numberOfParticles;
    int n = maxOrder;

    double t;
	double t_old;
    double delta_t;

    double absoluteError = absoluteErrorDefinition;
    double relativeError = relativeErrorDefinition;
    double infNorm;
    bool thisIsTrue;

    unsigned timeStepControlMethod=timeStepControlMethodDefinition;

    //double yPr1;
    //double yPr2;
    double distPr;
    double distPa;
	double r;
    double rmax;
    double rmin;

    double aMean=0.;
    double eMean=0.;

    double aM2=0.;
    double eM2=0.;

    double aVariance=0.;
    double eVariance=0.;

    unsigned rmaxCount=0;
    unsigned rminCount=0;

    double aDelta;
    double eDelta;

    double aNow;
    double eNow;

    unsigned periods=0;

    int pr=2;
    int pa=3;

    int mxi=rho_i;
    int mxj=rho_j;

    int coll=0;//collisions integer

    //double tolerance = regulaFalsiTolerance; //Convergence criterion tolerance

	double E0;
	double L0;

	double energyNow;
	double angMomNow;

    char buffer[15];

    std::ofstream File1;
    File1 << std::setprecision(16);
    sprintf(buffer,"P1_%.3u.dat",rank);
    File1.open(buffer);

    std::ofstream File2;
    File2 << std::setprecision(16);
    sprintf(buffer,"P2_%.3u.dat",rank);
    File2.open(buffer);

    std::ofstream File3;
    File3 << std::setprecision(16);
    sprintf(buffer,"P3_%.3u.dat",rank);
    File3.open(buffer);

    std::ofstream File4;
    File4 << std::setprecision(16);
    sprintf(buffer,"P4_%.3u.dat",rank);
    File4.open(buffer);

    std::ofstream File5;
    File5 << std::setprecision(16);
    sprintf(buffer,"P5_%.3u.dat",rank);
    File5.open(buffer);

    std::ofstream File6;
    File6 << std::setprecision(16);
    sprintf(buffer,"P6_%.3u.dat",rank);
    File6.open(buffer);

#ifdef FILE0
    std::ofstream File0;
    File0 << std::setprecision(16);
    sprintf(buffer,"P0_%.3u.dat",rank);
    File0.open(buffer);
#endif

    int tictac = 15*8*1986*29*rank+time(NULL);
    //int tictac = 1476959857;
	srand(tictac);
    std::cout << "SEED=" << tictac << std::endl;

    double rd1;
    double rd2;
    double rd3;
    double rd4;

    double a_RP;
    double e_RP;
    double theta_RP;
    double w_RP;

    double my_a;
    double my_e;

    for (unsigned c=1; c<=P; c++) {

        //std::cout << "c=" << c << " rank=" << rank << " ";

        blitz3Djet A_old (maxOrder,N+1,N+1);
        A_old=0.;
        //std::cout << "A_old(1,mxi,mxj)=" << A_old(1,mxi,mxj) << std::endl;

        blitz1Djet m(N+1);
        blitz2Djet x(maxOrder+1,N+1);
        blitz2Djet y(maxOrder+1,N+1);
        blitz2Djet x_ini(2,N+1);
        blitz2Djet y_ini(2,N+1);
        blitz3Djet X (maxOrder,N+1,N+1);
        blitz3Djet Y (maxOrder,N+1,N+1);
        blitz3Djet A (maxOrder,N+1,N+1); A=0.;
        blitz3Djet B (maxOrder,N+1,N+1);
        blitz3Djet C (maxOrder,N+1,N+1);
        blitz3Djet DX(maxOrder,N+1,N+1);
        blitz3Djet DY(maxOrder,N+1,N+1);
        blitz3Djet EX(maxOrder,N+1,N+1);
        blitz3Djet EY(maxOrder,N+1,N+1);

        m=M_Su,M_Sa,M_Ti,M_Pr,M_Pa,0.;

        t=0.;
        delta_t=timeStep;

        infNorm=0.;

        rd1 = ( 1.*rand() )/( 1.0*RAND_MAX );
        rd2 = ( 1.*rand() )/( 1.0*RAND_MAX );
        rd3 = ( 1.*rand() )/( 1.0*RAND_MAX );
        rd4 = ( 1.*rand() )/( 1.0*RAND_MAX );

//        //Condiciones iniciales que generan singularidad en met de pos fals,
//        //en variable m_...
//        a_RP = 2.014977953797247;
//        e_RP = 3.697961570554395e-05;
//        theta_RP = 3.380896879419916;
//        w_RP = 3.769077244345274;

        a_RP = a_RP_min+(a_RP_MAX-a_RP_min)*rd1;
        e_RP = e_RP_min+(e_RP_MAX-e_RP_min)*rd2;
        my_a=a_RP;
        my_e=e_RP;
        theta_RP = 360.*rd3;
        w_RP = 360.*rd4;

        newInitialConditions(x,y,m,index_sa,index_su,a_Sa,e_Sa,0.,0.);
        newInitialConditions(x,y,m,index_ti,index_sa,a_Ti,e_Ti,0.,0.);
        newInitialConditions(x,y,m,index_pr,index_sa,a_Pr,e_Pr,0.,0.);
        newInitialConditions(x,y,m,index_pa,index_sa,a_Pa,e_Pa,0.,0.);
        //InitialConditions(x,y,m,1,a_Ti,e_Ti,0.,0.);
        //InitialConditions(x,y,m,2,a_Pr,e_Pr,0.,0.);
        //InitialConditions(x,y,m,3,a_Pa,e_Pa,0.,0.);

        x(0,0)=(-m(1)*x(0,1)-m(2)*x(0,2)-m(3)*x(0,3)-m(4)*x(0,4))/m(0);
        y(0,0)=(-m(1)*y(0,1)-m(2)*y(0,2)-m(3)*y(0,3)-m(4)*y(0,4))/m(0);
        x(1,0)=(-m(1)*x(1,1)-m(2)*x(1,2)-m(3)*x(1,3)-m(4)*x(1,4))/m(0);
        y(1,0)=(-m(1)*y(1,1)-m(2)*y(1,2)-m(3)*y(1,3)-m(4)*y(1,4))/m(0);

        for (int i=0;i<N+1;i++){

            x_ini(0,i)=x(0,i);
            y_ini(0,i)=y(0,i);
            x_ini(1,i)=x(1,i);
            y_ini(1,i)=y(1,i);

        }

        //std::cout << "x00=" << x(0,0) << std::endl;
        //std::cout << "y00=" << y(0,0) << std::endl;
        //std::cout << "u00=" << x(1,0) << std::endl;
        //std::cout << "v00=" << y(1,0) << std::endl;

//        x(0,0)=0.;
//        y(0,0)=0.;
//        x(1,0)=0.;
//        y(1,0)=0.;

        E0 = Energy(N,x,y,m);
        L0 = AngularMomentum(N,x,y,m);
        if (c==1){

            std::cout << "Energy (blitz)=" << Energy(N,x,y,m) << std::endl;
            std::cout << "AngMom (blitz)=" << AngularMomentum(N,x,y,m) << std::endl;

        }

        std::cout << "c=" << c << " rank=" << rank << ", SEED=" << tictac << " ";

        newInitialConditions(x,y,m,index_rp,index_sa,a_RP,e_RP,theta_RP,w_RP);
        //InitialConditions(x,y,m,index_rp,a_RP,e_RP,theta_RP,w_RP);

        distPr=0.;
        distPa=0.;
        //yPr1=0.;
        //yPr2=0.;
        rmax=0.;
        rmin=0.;

        aMean=0.;
        eMean=0.;

        aM2=0.;
        eM2=0.;

        aVariance=0.;
        eVariance=0.;

        rmaxCount=0;
        rminCount=0;

        aDelta=0.;
        eDelta=0.;

        aNow=0.;
        eNow=0.;

        periods=0;

        coll=0;

		std::vector <double> MEANa;
		std::vector <double> MEANe;

        //for (unsigned d=0; periods<=maxperiods ; d++) {

        //if (c==282) maxperiods=2000;
        //else maxperiods=1;

        unsigned d;

        for (d=0; t/T_Pr<=maxperiods+1 ; d++) {

            //std::cout << "alright0" << std::endl;
            blitz2Djet xnew(maxOrder+1,N+1);
            blitz2Djet ynew(maxOrder+1,N+1);

            //std::cout << "alright01" << std::endl;
            xnew=x;
            ynew=y;
            //std::cout << "alright02" << std::endl;

            if (periods==maxperiods && coll==0) {

				//std::cout << "MEANa.size()=" << MEANa.size() << ", periods=" << periods << ", rminCount=" << rminCount << std::endl;
				double ss_a=0;
				double ss_e=0;
				for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
					ss_a=ss_a+MEANa[ss_i];
					ss_e=ss_e+MEANe[ss_i];
				}
				ss_a=ss_a/MEANa.size();
				ss_e=ss_e/MEANe.size();
				double vr_a=0;
				double vr_e=0;
				for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
					vr_a=vr_a+pow(MEANa[ss_i]-ss_a,2.);
					vr_e=vr_e+pow(MEANe[ss_i]-ss_e,2.);
				}
				vr_a=sqrt(vr_a/MEANa.size());
				vr_e=sqrt(vr_e/MEANe.size());

                File5 <<//
                        a_RP << " " <<//1
                        e_RP << " " <<//2
                        theta_RP << " " <<//3
                        w_RP << " " <<//4
                        ss_a << " " <<//5
                        ss_e << " " <<//6
                        vr_a << " " <<//7
                        vr_e << " " <<//8
                        t/T_Pr     << " " <<//9
                        xnew(0,0) << " " <<//10
                        ynew(0,0) << " " <<//11
                        xnew(1,0) << " " <<//12
                        ynew(1,0) << " " <<//13
                        xnew(0,1) << " " <<//14
                        ynew(0,1) << " " <<//15
                        xnew(1,1) << " " <<//16
                        ynew(1,1) << " " <<//17
                        xnew(0,2) << " " <<//18
                        ynew(0,2) << " " <<//19
                        xnew(1,2) << " " <<//20
                        ynew(1,2) << " " <<//21
                        xnew(0,3) << " " <<//22
                        ynew(0,3) << " " <<//23
                        xnew(1,3) << " " <<//24
                        ynew(1,3) << " " <<//25
                        xnew(0,4) << " " <<//26
                        ynew(0,4) << " " <<//27
                        xnew(1,4) << " " <<//28
                        ynew(1,4) << " " <<//29
                        xnew(0,5) << " " <<//30
                        ynew(0,5) << " " <<//31
                        xnew(1,5) << " " <<//32
                        ynew(1,5) << " " <<//33
                        my_a << " " <<//34
                        my_e << " " <<//35
                        aNow << " " <<//36
                        eNow << " " <<//37
                        x_ini(0,index_rp) << " " <<//38
                        y_ini(0,index_rp) << " " <<//39
                        x_ini(1,index_rp) << " " <<//40
                        y_ini(1,index_rp) << " " <<//41
                        (Energy(N,ynew,xnew,m)-E0)/E0 <<//42
                        std::endl;

                break;

                //std::cout << "periods=" << periods << std::endl;

            };//end, File5

            //yPr1=ynew(0,2);

            distPr=1.;
            distPa=1.;
            r=radialDistance(xnew,ynew,index_rp,index_sa);//0.5*(a_Pr+a_Pa);

            //std::cout << "alright" << std::endl;

//			double f_calc=myArcTan(x_rel,y_rel);
//			double v2_rel=u_rel*u_rel+v_rel*v_rel;
//			double r_rel=sqrt(x_rel*x_rel+y_rel*y_rel);
//			double E_calc=0.5*v2_rel-(G*m(0)/r_rel)*(1.+(0.5*J2)/(r_rel*r_rel));
//			double a_0=0.5*G/fabs(E_calc);
//			double a_calc=a_0-Gamma(E_calc,r_rel,a_0);
//			a_calc=a_calc-Gamma(E_calc,r_rel,a_calc);
//			a_calc=a_calc-Gamma(E_calc,r_rel,a_calc);
//			a_calc=a_calc-Gamma(E_calc,r_rel,a_calc);
//			a_calc=a_calc-Gamma(E_calc,r_rel,a_calc);
//			a_calc=a_calc-Gamma(E_calc,r_rel,a_calc);

#ifdef FILE0
            if(c==282){

				energyNow=Energy(N,x,y,m);
				angMomNow=AngularMomentum(N,x,y,m);

                File0 << t/T_Pr
                      << " " << aNow
                      << " " << eNow
                      << " " << semiMajorAxis(N,x,y,m,index_rp,index_sa)
                      << " " << eccentricity(N,x,y,m,index_rp,index_sa)
                      << " " << rminCount
                      << " " << rmaxCount
                      << " " << aMean
                      << " " << eMean
                      << " " << (energyNow-E0)/E0
                      << " " << (angMomNow-L0)/L0
                      << " " << x(0,0)
                      << " " << y(0,0)
                      << " " << x(0,1)
                      << " " << y(0,1)
                      << " " << x(0,2)
                      << " " << y(0,2)
                      << " " << x(0,3)
                      << " " << y(0,3)
                      << " " << x(0,4)
                      << " " << y(0,4)
                      << " " << x(0,5)
                      << " " << y(0,5)
                      //<< " " << m(0)*x(0,0)+m(1)*x(0,1)+m(2)*x(0,2)+m(3)*x(0,3)+m(4)*x(0,4)
                      //<< " " << m(0)*y(0,0)+m(1)*y(0,1)+m(2)*y(0,2)+m(3)*y(0,3)+m(4)*y(0,4)
                      //<< std::endl;
                      << " " << A(0,mxi,mxj)
                      << " " << A(1,mxi,mxj)
                      << " " << A(2,mxi,mxj)
                      << std::endl;

            }
#endif

			t_old=t;
            //int magic=1202;
            //if (c==393 && d==magic) std::cout << "***A_old=" << A_old(1,mxi,mxj) << std::endl;
            A_old=A;
            //if (c==393 && d==magic) std::cout << "**+A_old=" << A_old(1,mxi,mxj) << std::endl;
            TaylorMethod(n,N,x,y,m,X,Y,A,B,C,DX,DY,thisIsTrue,infNorm,absoluteError,relativeError,timeStepControlMethod,delta_t,xnew,ynew,t,maxperiods);
            if (A(1,mxi,mxj)*A_old(1,mxi,mxj)<0.){
                //double dtans=delta_t-A(1,mxi,mxj)*delta_t/(A(1,mxi,mxj)-A_old(1,mxi,mxj));
                ////double dtans=0.;
                //if (c==393 && d==magic) std::cout << "A=" << A(1,mxi,mxj) << ", A_old=" << A_old(1,mxi,mxj) << ", x=" << A_old(2,mxi,mxj) << ", c=" << c << ", d=" << d << std::endl;

                double a_=0.;
                double b_=delta_t;
                //double c_=0.5*delta_t; //bisection
                double c_=(A(1,mxi,mxj)*a_-A_old(1,mxi,mxj)*b_)/(A(1,mxi,mxj)-A_old(1,mxi,mxj)); //false position
                //double d_=-1.; //Illinois
                double f_=oneJetDotHornerSum(n,c_,A_old,mxi,mxj);

                for(int myl=0; fabs(0.5*(b_-a_))>newtonRaphsonTolerance && fabs(f_)>newtonRaphsonTolerance && myl<2;myl++){

                    if      (f_*oneJetDotHornerSum(n,a_,A_old,mxi,mxj)>0.) a_=c_;
                    else if (f_*oneJetDotHornerSum(n,b_,A_old,mxi,mxj)>0.) b_=c_;

                    //c_=0.5*(a_+b_); //bisection
                    //d_=c_; //Illinois
                    //if(d_==c_) c_=(oneJetDotHornerSum(n,b_,A_old,mxi,mxj)*a_-0.5*oneJetDotHornerSum(n,a_,A_old,mxi,mxj)*b_)/(oneJetDotHornerSum(n,b_,A_old,mxi,mxj)-0.5*oneJetDotHornerSum(n,a_,A_old,mxi,mxj)); //Illinois
                    c_=(oneJetDotHornerSum(n,b_,A_old,mxi,mxj)*a_-oneJetDotHornerSum(n,a_,A_old,mxi,mxj)*b_)/(oneJetDotHornerSum(n,b_,A_old,mxi,mxj)-oneJetDotHornerSum(n,a_,A_old,mxi,mxj)); //false position
                    f_=oneJetDotHornerSum(n,c_,A_old,mxi,mxj);

                    //std::cout << "mylbm=" << myl << std::endl;

                }
                double dtans=c_;
                ///*if(fabs(oneJetDotHornerSum(n,dtans,A_old,mxi,mxj))>newtonRaphsonTolerance)*/ std::cout << "A1_old=" << oneJetDotHornerSum(n,dtans,A_old,mxi,mxj) << ", c=" << c << ", d=" << d << ", rank=" << rank << std::endl;
                for(int myl=0;myl<4;myl++){
                    dtans=dtans-newtonQuotientHornerSum(n,dtans,A_old,mxi,mxj);
                    //std::cout << "mylnr=" << myl << std::endl;
                }

                ///*if(fabs(oneJetDotHornerSum(n,dtans,A_old,mxi,mxj))>newtonRaphsonTolerance)*/ std::cout << "A1_old=" << oneJetDotHornerSum(n,dtans,A_old,mxi,mxj) << ", c=" << c << ", d=" << d << ", rank=" << rank << std::endl;
                //if (c==393 && d==magic) std::cout << "t*/dt=" << dtans/delta_t << ", A0=" << oneJetHornerSum(n,dtans,A_old,mxi,mxj) << ", A1=" << oneJetDotHornerSum(n,dtans,A_old,mxi,mxj) << ", A2=" << oneJetDotDotHornerSum(n,dtans,A_old,mxi,mxj) << ", t/T_Pr=" << t/T_Pr << std::endl;
                //if (/*A_old(2,mxi,mxj)*/oneJetDotDotHornerSum(n,dtans,A_old,mxi,mxj)>0.){
                if (A_old(1,mxi,mxj)<0. && A(1,mxi,mxj)>0.){
                    rmin= sqrt(oneJetHornerSum(n,dtans,A_old,mxi,mxj));
                    //std::cout << "rmin=" << rmin << std::endl;
                    rminCount++;
                    //std::cout << "rmin=" << rmin << ", rminCount=" << rminCount << ", d=" << d << std::endl;
                    if( isnan(rmin) || isinf(rmin) ){

                        std::cout << "ISMIN rmin=" << rmin << ", c=" << c << ", d=" << d << ", rank=" << rank << ", SEED=" << tictac << std::endl;
                        std::cout << "t/T_Pr" << t/T_Pr << std::endl;
                        std::cout << "distPr_div_rHill=" << distPr/R_Hill_Pr << std::endl;
                        std::cout << "distPa_div_rHill=" << distPa/R_Hill_Pa << std::endl;
                        std::cout << "rmaxCount=" << rmaxCount << std::endl;
                        std::cout << "rminCount=" << rminCount << std::endl;
                        std::cout << "rminMean=" << aMean*(1.-eMean) << std::endl;
                        //rmin=aMean*(1.-eMean);

                    }

                }
                //else if (/*A_old(2,mxi,mxj)*/oneJetDotDotHornerSum(n,dtans,A_old,mxi,mxj)<0.){
                else if (A_old(1,mxi,mxj)>0. && A(1,mxi,mxj)<0.){
                    rmax= sqrt(oneJetHornerSum(n,dtans,A_old,mxi,mxj));
                    //std::cout << "rmax=" << rmax << std::endl;
                    rmaxCount++;
                    //std::cout << "rmax=" << rmax << ", rmaxCount=" << rmaxCount << ", d=" << d << std::endl;
                    if( isnan(rmax) || isinf(rmax) ){

                        std::cout << "ISMAX rmax=" << rmax << ", c=" << c << ", d=" << d << ", rank=" << rank << ", SEED=" << tictac << std::endl;
                        std::cout << "distPr_div_rHill=" << distPr/R_Hill_Pr << std::endl;
                        std::cout << "distPa_div_rHill=" << distPa/R_Hill_Pa << std::endl;
                        std::cout << "rmaxCount=" << rmaxCount << std::endl;
                        std::cout << "rminCount=" << rminCount << std::endl;
                        std::cout << "rminMean=" << aMean*(1.+eMean) << std::endl;
                        //rmax=aMean*(1.+eMean);

                    }

                }
                else std::cout << "WARNING: NR ill condition" << std::endl;

                /**/if(fabs((1.*rmaxCount)-(1.*rminCount))>1.)/**/ std::cout << "rMC=" << rmaxCount << ", rmC=" << rminCount << ", diff=" << fabs(1.*rmaxCount-1.*rminCount) << ", c=" << c << ", d=" << d << ", rank=" << rank << ", SEED=" << tictac << ", t/T=" << t/T_Pr << std::endl;

                if ( rmaxCount==rminCount ) {

                    //std::cout << "rmin=" << rmin << ", rmax=" << rmax << ", a=" << 0.5*(rmin+rmax) << ", e=" << (rmax-rmin)/(rmin+rmax) << std::endl;

                    aNow=0.5*(rmax+rmin);
                    if (rmax==rmin) eNow=0.;
                    else eNow=(rmax-rmin)/(2.*aNow);

                    if(rmaxCount==1 && rminCount==1 && periods<10 ){

                        a_RP=aNow;
                        e_RP=eNow;

                    }

                    aDelta = aNow - aMean;
                    aMean = aMean + aDelta/(1.*rmaxCount);
                    aM2 = aM2 + aDelta*(aNow-aMean);
                    aVariance = aM2/(1.*rmaxCount - 1.);

                    eDelta = eNow - eMean;
                    eMean = eMean + eDelta/(1.*rminCount);
                    eM2 = eM2 + eDelta*(eNow-eMean);
                    eVariance = eM2/(1.*rminCount - 1.);

//					std::cout << "rmaxCount=" << rmaxCount << " ,rminCount=" << rminCount << std::endl;
//					std::cout << "aNow=" << aNow << " ,eNow=" << eNow << std::endl;
//					std::cout << "aMean=" << aMean << " ,eMean=" << eMean << std::endl;
//					std::cout << "aM2=" << aM2 << " ,eM2=" << eM2 << std::endl;
//					std::cout << "aDelta=" << aDelta << " ,eDelta=" << eDelta << std::endl;
//					std::cout << "aVariance=" << aVariance << " ,eVariance=" << eVariance << std::endl;

                    if(rmaxCount==200 && rminCount==200){

                        MEANa.push_back(aMean);
                        MEANe.push_back(eMean);

                        rmaxCount=0;
                        rminCount=0;
                        aMean=0;
                        eMean=0;
                        aM2=0;
                        eM2=0;
                        eDelta=0;
                        aDelta=0;
                        aVariance=0;
                        eVariance=0;

                    }

                };//end, if rmaxCount==rminCount


            };//end, if A(1,mxi,mxj)*A_old(1,mxi,mxj)<0.

            //std::cout << "alright2" << std::endl;

            //yPr2=ynew(0,2);

            //std::cout << "alright3" << std::endl;

            //"Paint" particle if it enters Prometheus/Pandora Hill Sphere, or escapes the region between them:
            if ( distPr <= R_Hill_Pr || distPa <= R_Hill_Pa  || r>a_Pa || r<a_Pr || isnan(xnew(0,index_rp)) || isinf(xnew(0,index_rp))) {

                //std::cout << "         *rank = " << rank << " : Pro Coll!" << std::endl;
                //if (distPr <= R_Hill_Pr)        std::cout << "distPr/R_Hill_Pr="  << distPr/R_Hill_Pr  << std::endl;
                //else if (distPa <= R_Hill_Pa)   std::cout << "distPa/R_Hill_Pa="  << distPa/R_Hill_Pa << std::endl;
                //else                            std::cout << "*OUT OF BOUNDS* r=" << r << ", c=" << c << ", d=" << d << std::endl;

                coll=1;

                File6 <<//
                        a_RP << " " <<//
                        e_RP << " " <<//
                        theta_RP << " " <<//
                        w_RP << " " <<//
                        aMean << " " <<//
                        eMean << " " <<//
                        aVariance << " " <<//
                        eVariance << " " <<//
                        t/T_Pr     << " " <<//
                        xnew(0,0) << " " <<//
                        ynew(0,0) << " " <<//
                        xnew(1,0) << " " <<//
                        ynew(1,0) << " " <<//
                        xnew(0,1) << " " <<//
                        ynew(0,1) << " " <<//
                        xnew(1,1) << " " <<//
                        ynew(1,1) << " " <<//
                        xnew(0,2) << " " <<//
                        ynew(0,2) << " " <<//
                        xnew(1,2) << " " <<//
                        ynew(1,2) << " " <<//
                        xnew(0,3) << " " <<//
                        ynew(0,3) << " " <<//
                        xnew(1,3) << " " <<//
                        ynew(1,3) << " " <<//
                        xnew(0,4) << " " <<//
                        ynew(0,4) << " " <<//
                        xnew(1,4) << " " <<//
                        ynew(1,4) << " " <<//
                        xnew(0,5) << " " <<//
                        ynew(0,5) << " " <<//
                        xnew(1,5) << " " <<//
                        ynew(1,5) << " " <<//
                        my_a << " " <<//
                        my_e << " " <<//
                        aNow << " " <<//
                        eNow << " " <<//
                        x_ini(0,index_rp) << " " <<//38
                        y_ini(0,index_rp) << " " <<//39
                        x_ini(1,index_rp) << " " <<//40
                        y_ini(1,index_rp) << " " <<//41
                        (Energy(N,ynew,xnew,m)-E0)/E0 <<//
                        std::endl;

				break;

			};

            //if (yPr1<0 && yPr2>0) {
            if ( floor(t_old/T_Pr)<floor(t/T_Pr) ) {
                periods++;

                if (periods==snapshots1 && coll==0) {

					//std::cout << "MEANa.size()=" << MEANa.size() << ", periods=" << periods << ", rminCount=" << rminCount << std::endl;
					double ss_a=0;
					double ss_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						ss_a=ss_a+MEANa[ss_i];
						ss_e=ss_e+MEANe[ss_i];
					}
					ss_a=ss_a/MEANa.size();
					ss_e=ss_e/MEANe.size();
					double vr_a=0;
					double vr_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						vr_a=vr_a+pow(MEANa[ss_i]-ss_a,2.);
						vr_e=vr_e+pow(MEANe[ss_i]-ss_e,2.);
					}
					vr_a=sqrt(vr_a/MEANa.size());
					vr_e=sqrt(vr_e/MEANe.size());

                    File1 <<//
                            a_RP << " " <<//
                            e_RP << " " <<//
                            theta_RP << " " <<//
                            w_RP << " " <<//
                            ss_a << " " <<//
                            ss_e << " " <<//
                            vr_a << " " <<//
                            vr_e << " " <<//
                            t/T_Pr     << " " <<//
                            xnew(0,0) << " " <<//
                            ynew(0,0) << " " <<//
                            xnew(1,0) << " " <<//
                            ynew(1,0) << " " <<//
                            xnew(0,1) << " " <<//
                            ynew(0,1) << " " <<//
                            xnew(1,1) << " " <<//
                            ynew(1,1) << " " <<//
                            xnew(0,2) << " " <<//
                            ynew(0,2) << " " <<//
                            xnew(1,2) << " " <<//
                            ynew(1,2) << " " <<//
                            xnew(0,3) << " " <<//
                            ynew(0,3) << " " <<//
                            xnew(1,3) << " " <<//
                            ynew(1,3) << " " <<//
                            xnew(0,4) << " " <<//
                            ynew(0,4) << " " <<//
                            xnew(1,4) << " " <<//
                            ynew(1,4) << " " <<//
                            xnew(0,5) << " " <<//
                            ynew(0,5) << " " <<//
                            xnew(1,5) << " " <<//
                            ynew(1,5) << " " <<//
                            my_a << " " <<//
                            my_e << " " <<//
                            aNow << " " <<//
                            eNow << " " <<//
                            x_ini(0,index_rp) << " " <<//38
                            y_ini(0,index_rp) << " " <<//39
                            x_ini(1,index_rp) << " " <<//40
                            y_ini(1,index_rp) << " " <<//41
                            (Energy(N,ynew,xnew,m)-E0)/E0 <<//
                            std::endl;

                }//end, File1
                else if (periods==snapshots2 && coll==0) {

					//std::cout << "MEANa.size()=" << MEANa.size() << ", periods=" << periods << ", rminCount=" << rminCount << std::endl;
					double ss_a=0;
					double ss_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						ss_a=ss_a+MEANa[ss_i];
						ss_e=ss_e+MEANe[ss_i];
					}
					ss_a=ss_a/MEANa.size();
					ss_e=ss_e/MEANe.size();
					double vr_a=0;
					double vr_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						vr_a=vr_a+pow(MEANa[ss_i]-ss_a,2.);
						vr_e=vr_e+pow(MEANe[ss_i]-ss_e,2.);
					}
					vr_a=sqrt(vr_a/MEANa.size());
					vr_e=sqrt(vr_e/MEANe.size());

                    File2 <<//
                            a_RP << " " <<//
                            e_RP << " " <<//
                            theta_RP << " " <<//
                            w_RP << " " <<//
                            ss_a << " " <<//
                            ss_e << " " <<//
                            vr_a << " " <<//
                            vr_e << " " <<//
                            t/T_Pr     << " " <<//
                            xnew(0,0) << " " <<//
                            ynew(0,0) << " " <<//
                            xnew(1,0) << " " <<//
                            ynew(1,0) << " " <<//
                            xnew(0,1) << " " <<//
                            ynew(0,1) << " " <<//
                            xnew(1,1) << " " <<//
                            ynew(1,1) << " " <<//
                            xnew(0,2) << " " <<//
                            ynew(0,2) << " " <<//
                            xnew(1,2) << " " <<//
                            ynew(1,2) << " " <<//
                            xnew(0,3) << " " <<//
                            ynew(0,3) << " " <<//
                            xnew(1,3) << " " <<//
                            ynew(1,3) << " " <<//
                            xnew(0,4) << " " <<//
                            ynew(0,4) << " " <<//
                            xnew(1,4) << " " <<//
                            ynew(1,4) << " " <<//
                            xnew(0,5) << " " <<//
                            ynew(0,5) << " " <<//
                            xnew(1,5) << " " <<//
                            ynew(1,5) << " " <<//
                            my_a << " " <<//
                            my_e << " " <<//
                            aNow << " " <<//
                            eNow << " " <<//
                            x_ini(0,index_rp) << " " <<//38
                            y_ini(0,index_rp) << " " <<//39
                            x_ini(1,index_rp) << " " <<//40
                            y_ini(1,index_rp) << " " <<//41
                            (Energy(N,ynew,xnew,m)-E0)/E0 <<//
                            std::endl;

                    //std::cout << "periods=" << periods << std::endl;

                }//end, File2
                else if (periods==snapshots3 && coll==0) {

					//std::cout << "MEANa.size()=" << MEANa.size() << ", periods=" << periods << ", rminCount=" << rminCount << std::endl;
					double ss_a=0;
					double ss_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						ss_a=ss_a+MEANa[ss_i];
						ss_e=ss_e+MEANe[ss_i];
					}
					ss_a=ss_a/MEANa.size();
					ss_e=ss_e/MEANe.size();
					double vr_a=0;
					double vr_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						vr_a=vr_a+pow(MEANa[ss_i]-ss_a,2.);
						vr_e=vr_e+pow(MEANe[ss_i]-ss_e,2.);
					}
					vr_a=sqrt(vr_a/MEANa.size());
					vr_e=sqrt(vr_e/MEANe.size());

                    File3 <<//
                            a_RP << " " <<//
                            e_RP << " " <<//
                            theta_RP << " " <<//
                            w_RP << " " <<//
                            ss_a << " " <<//
                            ss_e << " " <<//
                            vr_a << " " <<//
                            vr_e << " " <<//
                            t/T_Pr     << " " <<//
                            xnew(0,0) << " " <<//
                            ynew(0,0) << " " <<//
                            xnew(1,0) << " " <<//
                            ynew(1,0) << " " <<//
                            xnew(0,1) << " " <<//
                            ynew(0,1) << " " <<//
                            xnew(1,1) << " " <<//
                            ynew(1,1) << " " <<//
                            xnew(0,2) << " " <<//
                            ynew(0,2) << " " <<//
                            xnew(1,2) << " " <<//
                            ynew(1,2) << " " <<//
                            xnew(0,3) << " " <<//
                            ynew(0,3) << " " <<//
                            xnew(1,3) << " " <<//
                            ynew(1,3) << " " <<//
                            xnew(0,4) << " " <<//
                            ynew(0,4) << " " <<//
                            xnew(1,4) << " " <<//
                            ynew(1,4) << " " <<//
                            xnew(0,5) << " " <<//
                            ynew(0,5) << " " <<//
                            xnew(1,5) << " " <<//
                            ynew(1,5) << " " <<//
                            my_a << " " <<//
                            my_e << " " <<//
                            aNow << " " <<//
                            eNow << " " <<//
                            x_ini(0,index_rp) << " " <<//38
                            y_ini(0,index_rp) << " " <<//39
                            x_ini(1,index_rp) << " " <<//40
                            y_ini(1,index_rp) << " " <<//41
                            (Energy(N,ynew,xnew,m)-E0)/E0 <<//
                            std::endl;

                    //std::cout << "periods=" << periods << std::endl;

                }//end, File3
                else if (periods==snapshots4 && coll==0) {

					//std::cout << "MEANa.size()=" << MEANa.size() << ", periods=" << periods << ", rminCount=" << rminCount << std::endl;
					double ss_a=0;
					double ss_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						ss_a=ss_a+MEANa[ss_i];
						ss_e=ss_e+MEANe[ss_i];
					}
					ss_a=ss_a/MEANa.size();
					ss_e=ss_e/MEANe.size();
					double vr_a=0;
					double vr_e=0;
					for (int ss_i=0; ss_i<MEANa.size(); ss_i++) {
						vr_a=vr_a+pow(MEANa[ss_i]-ss_a,2.);
						vr_e=vr_e+pow(MEANe[ss_i]-ss_e,2.);
					}
					vr_a=sqrt(vr_a/MEANa.size());
					vr_e=sqrt(vr_e/MEANe.size());

                    File4 <<//
                            a_RP << " " <<//
                            e_RP << " " <<//
                            theta_RP << " " <<//
                            w_RP << " " <<//
                            ss_a << " " <<//
                            ss_e << " " <<//
                            vr_a << " " <<//
                            vr_e << " " <<//
                            t/T_Pr     << " " <<//
                            xnew(0,0) << " " <<//
                            ynew(0,0) << " " <<//
                            xnew(1,0) << " " <<//
                            ynew(1,0) << " " <<//
                            xnew(0,1) << " " <<//
                            ynew(0,1) << " " <<//
                            xnew(1,1) << " " <<//
                            ynew(1,1) << " " <<//
                            xnew(0,2) << " " <<//
                            ynew(0,2) << " " <<//
                            xnew(1,2) << " " <<//
                            ynew(1,2) << " " <<//
                            xnew(0,3) << " " <<//
                            ynew(0,3) << " " <<//
                            xnew(1,3) << " " <<//
                            ynew(1,3) << " " <<//
                            xnew(0,4) << " " <<//
                            ynew(0,4) << " " <<//
                            xnew(1,4) << " " <<//
                            ynew(1,4) << " " <<//
                            xnew(0,5) << " " <<//
                            ynew(0,5) << " " <<//
                            xnew(1,5) << " " <<//
                            ynew(1,5) << " " <<//
                            my_a << " " <<//
                            my_e << " " <<//
                            aNow << " " <<//
                            eNow << " " <<//
                            x_ini(0,index_rp) << " " <<//38
                            y_ini(0,index_rp) << " " <<//39
                            x_ini(1,index_rp) << " " <<//40
                            y_ini(1,index_rp) << " " <<//41
                            (Energy(N,ynew,xnew,m)-E0)/E0 <<//
                            std::endl;

                    //std::cout << "periods=" << periods << std::endl;

                };//end, File4

                //std::cout << "periods=" << periods << std::endl;

            }; //end, if floor(t_old/T_Pr)<floor(t/T_Pr)

            //std::cout << "alright4" << std::endl;
            x=xnew;
            y=ynew;
            //std::cout << "alright5" << std::endl;

            //

        };//end, for d

        std::cout << "t/T_Pr=" << t/T_Pr << ", d=" << d << std::endl;
        //if(c==numberOfInitialConditions)std::cout << "c=" << c << std::endl;

#ifdef FILE0
            File0 << std::endl;
#endif

    };//end, for c

    std::cout << "Done." << std::endl;

#ifdef PARALLEL
    std::cout << " Process #" << rank << " on " << processor_name << "." << std::endl;
	MPI_Finalize();
#endif

#ifdef FILE0
	File0.close();
#endif

    File1.close();
    File2.close();
    File3.close();
    File4.close();
    File5.close();
    File6.close();

    return 0;

}
