//
//  blitztest.h
//
//
//  Created by Jorge Antonio Perez Hernandez on 5/13/14.
//
//

#ifndef ____blitztest__
#define ____blitztest__

#include <iostream>

#endif /* defined(____blitztest__) */

//double Gamma(double E, double r, double a){
//
//	double A=(fabs(E)/G-Lambda/(r*r*r))*2.;
//	double C=6.*Lambda/r;
//	double D=-3.*Lambda;
//
//	return (((A*a-1.)*a+C)*a+D)/((3.*A*a-2.)*a+C);
//
//}
//
//double myArcTan(double x,double y){
//
//	if (x>0) return atan(y/x);
//	else if (x<0){
//		if (y>0) return atan(y/x)+pi;
//		else if (y<0) return atan(y/x)-pi;
//		else return -1./0.;
//	}
//	else return 1./0.;
//}

//double delta(int &N,
//             blitz2Djet &x,
//             blitz2Djet &y,
//             blitz1Djet &m,
//             int i,
//             int sign
//             ){

//    double my_hz=(x(0,i)-x(0,0))*(y(1,i)-y(1,0))-(y(0,i)-y(0,0))*(x(1,i)-x(1,0));

//    double my_mu=G*(1.+m(i));

//    double my_E=(my_hz*my_hz)/(2.*())

//    return pow(  (2./relativeDistance)-(relativeSpeed*relativeSpeed)/(G*(1.+m(i)))  ,-1.);

//}//end, delta

void newInitialConditions(blitz2Djet &x, //x,px
                          blitz2Djet &y, //y,py
                          blitz1Djet &m, //mass
                          int i, //ith particle
                          int prind, //primary index
                          double a, //semimajor axis
                          double e, //eccentricity
                          double ta, //true anomaly (Murray's f)
                          double lp //longitude of periapse
){

    double f=pi*ta/180.;
    double w=pi*lp/180.;

    double myMeanMotion = sqrt( G*(m(prind)+m(i))/(a*a*a) )*(1.+0.5*Lambda/(a*a)-0.125*Lambda*Lambda/(a*a*a*a)+2.*Lambda*e*e/(a*a));
    double myRadius = a*(1.-e*e)/(1.+e*cos(f));

    double x_0 = myRadius*cos(f);
    double y_0 = myRadius*sin(f);

    double u_0 = -myMeanMotion*a*   sin(f) /sqrt(1.-e*e);
    double v_0 =  myMeanMotion*a*(e+cos(f))/sqrt(1.-e*e);

    double nx0 = x_0*cos(w)-y_0*sin(w);
    double ny0 = x_0*sin(w)+y_0*cos(w);
    double nu0 = u_0*cos(w)-v_0*sin(w);
    double nv0 = u_0*sin(w)+v_0*cos(w);

    x(0,i)=nx0;
    y(0,i)=ny0;
    x(1,i)=nu0;
    y(1,i)=nv0;

//    std::cout << "x0"<< i <<"=" << nx0 << std::endl;
//    std::cout << "y0"<< i <<"=" << ny0 << std::endl;
//    std::cout << "u0"<< i <<"=" << nu0 << std::endl;
//    std::cout << "v0"<< i <<"=" << nv0 << std::endl;

    if (prind>0){
        x(0,i)+=x(0,prind);
        y(0,i)+=y(0,prind);
        x(1,i)+=x(1,prind);
        y(1,i)+=y(1,prind);
    }

    if(i==index_rp){

//        x(0,i)=22018.97874320764;
//        y(0,i)=4451.287903277223;
//        x(1,i)=-1.193883385786046;
//        y(1,i)=-0.7110904320632346;

        std::cout << "x0"<< i <<"=" << x(0,i) << " ";
//        std::cout << "y0"<< i <<"=" << y(0,i) << " ";
//        std::cout << "u0"<< i <<"=" << x(1,i) << " ";
//        std::cout << "v0"<< i <<"=" << y(1,i) << " ";

    }


}

double oneJetHornerSum(int &n,
                       double my_delta_t,
                       blitz3Djet &jet,
                       int my_i,
                       int my_j
                       ){

    double sum=jet(n,my_i,my_j);

    for(int l=n-1;l>=0;l--){
        sum=jet(l,my_i,my_j)+sum*my_delta_t;
    }

    return sum;

}

double oneJetDotHornerSum(int &n,
                       double my_delta_t,
                       blitz3Djet &jet,
                       int my_i,
                       int my_j
                       ){

    double sum=n*jet(n,my_i,my_j);

    for(int l=n-1;l>=1;l--){
        sum=l*jet(l,my_i,my_j)+sum*my_delta_t;
    }

    return sum;

}

double oneJetDotDotHornerSum(int &n,
                       double my_delta_t,
                       blitz3Djet &jet,
                       int my_i,
                       int my_j
                       ){

    double sum=n*(n-1)*jet(n,my_i,my_j);

    for(int l=n-1;l>=2;l--){
        sum=l*(l-1)*jet(l,my_i,my_j)+sum*my_delta_t;
    }

    return sum;

}

double newtonQuotientHornerSum(int &n,
                       double my_delta_t,
                       blitz3Djet &jet,
                       int my_i,
                       int my_j
                       ){

    double num=n*jet(n,my_i,my_j);
    double den=(n-1)*num;

    for(int l=n-1;l>=2;l--){
        num=l      *jet(l,my_i,my_j)+num*my_delta_t;
        den=l*(l-1)*jet(l,my_i,my_j)+den*my_delta_t;
    }

    num=1.*jet(1,my_i,my_j)+num*my_delta_t;

    return num/den;

}

double Energy(int &N,
              blitz2Djet &x,
              blitz2Djet &y,
              blitz1Djet &m
              ){

    double ans = 0.;

    for (int i=0; i<N; i++) {
        ans=ans+0.5*m(i)*(x(1,i)*x(1,i)+y(1,i)*y(1,i));

        for (int j=0; j<i; j++) {

            double rij=sqrt(  ( x(0,i)-x(0,j) )*( x(0,i)-x(0,j) )+( y(0,i)-y(0,j) )*( y(0,i)-y(0,j) )  );

            if (i==index_sa || j==index_sa) {
                ans=ans-G*m(j)*m(i)*(1.+Lambda/(3.*rij*rij))/rij;
            }
            else if(i!=index_sa && j!=index_sa) {
                ans=ans-G*m(j)*m(i)/rij;
            }
            else std::cout << "ERROR: Energy ERROR 1" << std::endl;

        };//end, for j

    };//end, for i

    return ans;

}//end, Energy

double AngularMomentum(int &N,
					   blitz2Djet &x,
					   blitz2Djet &y,
					   blitz1Djet &m
					   ){

	double ans = 0.;

	for (int i=0; i<N; i++) {
		ans=ans+m(i)*(x(0,i)*y(1,i)-y(0,i)*x(1,i));

	};//end, for i

	return ans;

}//end, AngularMomentum

double radialDistance(blitz2Djet &x,
                     blitz2Djet &y,
                     int i,
                     int j
                     ){

    return sqrt(  pow(x(0,i)-x(0,j),2.)+pow(y(0,i)-y(0,j),2.)  );

}//end, radialDistance (relative)

double radialVelocityNumerator(blitz2Djet &x,
                     blitz2Djet &y,
                     int i,
                     int j
                     ){

    return (x(0,i)-x(0,j))*(x(1,i)-x(1,j))+(y(0,i)-y(0,j))*(y(1,i)-y(1,j));

}//end, radialVelocityNumerator

double radialDistanceDotDot(blitz2Djet &x,
                     blitz2Djet &y,
                     int i,
                     int j
                     ){

    double term1= pow(x(1,i)-x(1,j),2.)+pow(y(1,i)-y(1,j),2.);
    double term2= (x(0,i)-x(0,j))*(x(2,i)-x(2,j))+(y(0,i)-y(0,j))*(y(2,i)-y(2,j));
    double term3= -pow(radialVelocityNumerator(x,y,i,j)/radialDistance(x,y,i,j),2.);

    return (term1+2.*term2+term3)/radialDistance(x,y,i,j);

}//end, radialDistanceDotDot

double unit2BP_hz(int &N,
                     blitz2Djet &x,
                     blitz2Djet &y,
                     blitz1Djet &m,
                     int i,
                     int j
                     ){

    return (x(0,i)-x(0,j))*(y(1,i)-y(1,j))-(y(0,i)-y(0,j))*(x(1,i)-x(1,j));

}//end, unit2BP_hz, Murray pg. 24 Eq. 2.6; pg. 25 Eq. 2.8

double LRLx(int &N,
			blitz2Djet &x,
			blitz2Djet &y,
			blitz1Djet &m,
            int i,
            int j
			){

    double relativeAngular=unit2BP_hz(N,x,y,m,i,j);

    double relativeDistance=sqrt(  pow((x(0,i)-x(0,j)),2.)  +  pow((y(0,i)-y(0,j)),2.)  );

    return +(y(1,i)-y(1,j))*relativeAngular/(G*(1.+m(i)))-(x(0,i)-x(0,j))/relativeDistance;

}//end, Laplace-Runge-Lenz vector, x component

double LRLy(int &N,
			blitz2Djet &x,
			blitz2Djet &y,
			blitz1Djet &m,
            int i,
            int j
			){

    double relativeAngular=unit2BP_hz(N,x,y,m,i,j);

    double relativeDistance=sqrt(  pow((x(0,i)-x(0,j)),2.)  +  pow((y(0,i)-y(0,j)),2.)  );

    return -(x(1,i)-x(1,j))*relativeAngular/(G*(1.+m(i)))-(y(0,i)-y(0,j))/relativeDistance;

}//end, Laplace-Runge-Lenz vector, y component

double semiMajorAxis(int &N,
					 blitz2Djet &x,
					 blitz2Djet &y,
					 blitz1Djet &m,
                     int i,
                     int j
					 ){

    double relativeDistance=sqrt(  pow((x(0,i)-x(0,j)),2.)  +  pow((y(0,i)-y(0,j)),2.)  );

    double relativeSpeed   =pow((x(1,i)-x(1,j)),2.)  +  pow((y(1,i)-y(1,j)),2.);

    return pow(  (2./relativeDistance)-(relativeSpeed)/(G*(1.+m(i)))  ,-1.);

}//end, semiMajorAxis, Murray pg. 53 Ec. 2.134

double eccentricity(int &N,
					blitz2Djet &x,
					blitz2Djet &y,
					blitz1Djet &m,
                    int i,
                    int j
					){

    double current_sma=semiMajorAxis(N,x,y,m,i,j);

    double current_ams=unit2BP_hz(N,x,y,m,i,j); //angular momentum per unit mass
	current_ams=current_ams*current_ams; //square it!

	return sqrt(  1.-current_ams/(G*(1.+m(i))*current_sma)  );

}//end, eccentricity, Murray pg. 53 Ec. 2.135

double sin_f(int &N,
			 blitz2Djet &x,
			 blitz2Djet &y,
			 blitz1Djet &m,
             int i,
             int j
			 ){

    double current_sma=semiMajorAxis(N,x,y,m,i,j);
    double current_ecc=eccentricity (N,x,y,m,i,j);
    double current_uhz=unit2BP_hz   (N,x,y,m,i,j);

    double relativeDistance=sqrt(  pow((x(0,i)-x(0,j)),2.)  +  pow((y(0,i)-y(0,j)),2.)  );

    double radialSpeed=(x(0,i)-x(0,j))*(x(1,i)-x(1,j))+(y(0,i)-y(0,j))*(y(1,i)-y(1,j));
	radialSpeed=radialSpeed/relativeDistance;

	return (  (current_sma*(1-current_ecc*current_ecc))/(current_uhz*current_ecc)  )*radialSpeed;

}//end, sin_f, Murray pg. 53 Ec. 2.139

double cos_f(int &N,
			 blitz2Djet &x,
			 blitz2Djet &y,
			 blitz1Djet &m,
             int i,
             int j
			 ){

    double current_sma=semiMajorAxis(N,x,y,m,i,j);
    double current_ecc=eccentricity (N,x,y,m,i,j);

    double relativeDistance=sqrt(  pow((x(0,i)-x(0,j)),2.)  +  pow((y(0,i)-y(0,j)),2.)  );

	return (1./current_ecc)*(  current_sma*(1.-current_ecc*current_ecc)/relativeDistance  -  1  );

}//end, cos_f, Murray pg. 53 Ec. 2.139

void InitialConditions(blitz2Djet &x, //x,px
                          blitz2Djet &y, //y,py
                          blitz1Djet &m, //mass
                          int i, //ith particle
                          double a, //semimajor axis
                          double e, //eccentricity
                          double theta, //true longitude
                          double w //longitude of periapse
){

    double n = sqrt( G*(1.+m(i))*( 1.+Lambda/(a*a) )/( a*a*a )  );
    double kappa = sqrt( G*(1.+m(i))*( 1.-Lambda/(a*a) )/( a*a*a )  );
//	double n = sqrt( G*( 1.+Lambda/(a*a) )/( a*a*a )  );

    //std::cout << "InitialConditions:: n=" << n << std::endl;

	double ai=a; //a/(1.+m(i));

    double r = ai*(1.0-e*e)/(1.0+e*cos(theta-w));

    double x0 = r*cos(theta);
    double y0 = r*sin(theta);

    double x_dot0 = -( n*ai*(e*sin(w)+sin(theta)) )/sqrt(1-e*e);
    double y_dot0 =  ( n*ai*(e*cos(w)+cos(theta)) )/sqrt(1-e*e);

    x(0,i)=x0;
    y(0,i)=y0;
    x(1,i)=x_dot0;
    y(1,i)=y_dot0;

//    std::cout << "INITIAL CONDITIONS" << std::endl;
    std::cout << "x0     = " << x0     << std::endl;
//    std::cout << "y0     = " << y0     << std::endl;
//    std::cout << "x_dot0 = " << x_dot0 << std::endl;
//    std::cout << "y_dot0 = " << y_dot0 << std::endl;

//    std::cout << "n      = " << n      << std::endl;
//    std::cout << "kappa  = " << kappa  << std::endl;
//    std::cout << "a      = " << a      << std::endl;
//    std::cout << "e      = " << e      << std::endl;

//    std::cout << "theta  = " << theta  << std::endl;
//    std::cout << "w      = " << w      << std::endl;

}//end, InitialConditions

void OrderControl(int &n,
                  int &N,
                  blitz2Djet &x,
                  blitz2Djet &y,
                  bool &thisIsTrue,
                  double &infNorm,
                  double &absoluteError,
                  double &relativeError,
                  unsigned &timeStepControlMethod
                  ){

    if (timeStepControlMethod == 1 || timeStepControlMethod == 2) {

        std::vector <double> allxy;

        //Building a vector which contains all position coordinates at present time (in abs value):
        for (int i=0;i<N+1;i++){
            allxy.push_back(  fabs( x(0,i) )  );
            allxy.push_back(  fabs( y(0,i) )  );
            //std::cout << "fabs( x[0]["<<i<<"] )=" << fabs( x[0][i] ) << std::endl;
            //std::cout << "fabs( y[0]["<<i<<"] )=" << fabs( y[0][i] ) << std::endl;
        };//end, for i

        //Calculating the absolutely-largest element of this vector:
        std::vector<double>::iterator normInf;
        normInf = max_element( allxy.begin() , allxy.end() );

        infNorm=*normInf;

        if (relativeError*infNorm <= absoluteError) {
            n = ceil( 1.0-0.5*log(absoluteError) );
            thisIsTrue = true;
        }//end, if relativeError*infNorm <= absoluteError
        else {
            n = ceil( 1.0-0.5*log(relativeError) );
            thisIsTrue = false;
        };//end, if relativeError*infNorm > absoluteError

        if (n>40) n=40;

    };//end, if timeStepControlMethod == 1 || timeStepControlMethod == 2

    //std::cout << "n=" << n << std::endl;

}//end, OrderControl

void TimeStepControl(int &n,
                     int &N,
                     blitz2Djet &x,
                     blitz2Djet &y,
                     bool &thisIsTrue,
                     double &infNorm,
                     double &absoluteError,
                     double &relativeError,
                     unsigned &timeStepControlMethod,
                     double &delta_t
                     ){

    if (timeStepControlMethod == 1 || timeStepControlMethod == 2 ) {

        //Calculating optimal time step:
        double rho_n;
        double rho_nm1;

        double zeta; //2nd step size control!

        std::vector <double> n_th_Jet; //   n   th degree Taylor Expansion Coefficients,
        std::vector <double> nm1_th_Jet; // n-1 th degree Taylor Expansion Coefficients.

        for (int i=0;i<N+1;i++){
            n_th_Jet  .push_back( fabs(x(n  ,i)) );
            n_th_Jet  .push_back( fabs(y(n  ,i)) );
            nm1_th_Jet.push_back( fabs(x(n-1,i)) );
            nm1_th_Jet.push_back( fabs(y(n-1,i)) );
        }//end, for i

        std::vector <double>::iterator it_for_n__;
        std::vector <double>::iterator it_for_nm1;

        //Calculating the absolutely-largest element of nth and n-1th jet:
        it_for_n__ = max_element( n_th_Jet  .begin() , n_th_Jet  .end() );
        it_for_nm1 = max_element( nm1_th_Jet.begin() , nm1_th_Jet.end() );

        double infty_norm_n_th_Jet   = *it_for_n__;
        double infty_norm_nm1_th_Jet = *it_for_nm1;

        if (thisIsTrue) {
            rho_n  = pow(  infty_norm_n_th_Jet    ,  ( (-1.0)/(1.0* n   ) )  );
            rho_nm1= pow(  infty_norm_nm1_th_Jet  ,  ( (-1.0)/(1.0*(n-1)) )  );
            zeta = 1.;
        }
        else {
            rho_n  = pow(  infNorm/infty_norm_n_th_Jet    ,  ( (1.0)/(1.0* n   ) )  );
            rho_nm1= pow(  infNorm/infty_norm_nm1_th_Jet  ,  ( (1.0)/(1.0*(n-1)) )  );
            zeta = infNorm;
        };//end, if thisIsTrue == true :P

        //1st step size control method: safety factor
        delta_t = fmin( rho_n , rho_nm1 )*exp(-2.)*exp(  (safetyFactor/( 1.0*(n-1) ))  );

        //2nd step size control method: succesive series terms cannot increase!
        if (timeStepControlMethod == 2){

            blitz1Djet infNormsAllOrders(n+1);
            unsigned zetaBadCount=0;

            for(int k=0;k<=n;k++){

                std::vector <double> k_th_Jet;

                for (int j=0;j<N+1;j++){
                    k_th_Jet.push_back( fabs(x(k,j)) );
                    k_th_Jet.push_back( fabs(y(k,j)) );
                };//end, for j

                std::vector <double>::iterator it_for_k;
                it_for_k = max_element( k_th_Jet.begin() , k_th_Jet.end() );

                infNormsAllOrders(k) = *it_for_k;

                if ( infNormsAllOrders(k)*pow(delta_t,1.*k) > zeta )zetaBadCount++;

            };//end, for k

            double newdelta_t=delta_t;

            for (unsigned l=1; zetaBadCount!=0 && l<20; l++) {
                //std::cout << "                    ********l=" << l << std::endl;
                zetaBadCount=0;
                for (unsigned k=0; k<=n; k++) {
                    if ( infNormsAllOrders(k)*pow(newdelta_t,1.*k) > zeta )zetaBadCount++;
                };//end, for k
                if (zetaBadCount>0) newdelta_t=exp(-1.*l)*delta_t;
            };//end, for l

            delta_t=newdelta_t;

        };//end, if method==2

    };//end, if timeStepControlMethod == 1 || timeStepControlMethod == 2

    //std::cout << "delta_t/tpr=" << delta_t/T_Pr << std::endl;

}//end, TimeStepControl

void TaylorSeriesExpansion(int &n,
                           int &N,
                           blitz2Djet &x,
                           blitz2Djet &y,
                           blitz1Djet &m,
                           double &delta_t,
                           blitz2Djet &xnew,
                           blitz2Djet &ynew
                           ){

    for (int j=1; j<N+1; j++) {

        double xnew0j=  x(n,j);
        double ynew0j=  y(n,j);
        double xnew1j=n*x(n,j);
        double ynew1j=n*y(n,j);

        for (int k=n-1; k>=1; k--) {
            xnew0j=xnew0j*delta_t+  x(k,j);
            ynew0j=ynew0j*delta_t+  y(k,j);
            xnew1j=xnew1j*delta_t+k*x(k,j);
            ynew1j=ynew1j*delta_t+k*y(k,j);
        };//end, for k

        xnew(0,j)=(xnew0j*delta_t+x(0,j));
        ynew(0,j)=(ynew0j*delta_t+y(0,j));
        xnew(1,j)=xnew1j;
        ynew(1,j)=ynew1j;

    };//end, for j

    xnew(0,0)=( -m(1)*xnew(0,1)-m(2)*xnew(0,2)-m(3)*xnew(0,3) )/m(0);
    ynew(0,0)=( -m(1)*ynew(0,1)-m(2)*ynew(0,2)-m(3)*ynew(0,3) )/m(0);
    xnew(1,0)=( -m(1)*xnew(1,1)-m(2)*xnew(1,2)-m(3)*xnew(1,3) )/m(0);
    ynew(1,0)=( -m(1)*ynew(1,1)-m(2)*ynew(1,2)-m(3)*ynew(1,3) )/m(0);

    //std::cout << "n=" << n << std::endl;
    //std::cout << "delta_t/tpr=" << delta_t/T_Pr << std::endl;

}//end, TaylorSeriesExpansion

void TaylorMethod(int &n,
                  int &N,
                  blitz2Djet &x,
                  blitz2Djet &y,
                  blitz1Djet &m,
                  blitz3Djet &X,
                  blitz3Djet &Y,
                  blitz3Djet &A,
                  blitz3Djet &B,
                  blitz3Djet &C,
                  blitz3Djet &DX,
                  blitz3Djet &DY,
                  bool &thisIsTrue,
                  double &infNorm,
                  double &absoluteError,
                  double &relativeError,
                  unsigned &timeStepControlMethod,
                  double &delta_t,
                  blitz2Djet &xnew,
                  blitz2Djet &ynew,
                  double &t,
                  int maxperiods
                  ){

    OrderControl(n,N,x,y,thisIsTrue,infNorm,absoluteError,relativeError,timeStepControlMethod);

    double oddOrder,evenOrder,xorderi,yorderi;
    int parity;

    //std::cout << "SEARCHING" << std::endl;
    for (int order=0; order<=n-2; order++) {

        oddOrder = (order-1.)/2.;
        evenOrder = (order/2.)-1.;
        parity    = order%2;

        //        std::cout << "    order=" <<     order << std::endl;
        //        std::cout << " oddOrder=" <<  oddOrder << std::endl;
        //        std::cout << "evenOrder=" << evenOrder << std::endl;
        //        std::cout << "   parity=" <<    parity << std::endl;

        double ans,bns,cns,dns;

        for (int i=0 ; i<N ; ++i){

            xorderi=x(order,i);
            yorderi=y(order,i);

            for (int j = i+1; j < N+1; ++ j){

                //X,Y-Jet calculation:
                X(order,i,j)=xorderi-x(order,j);
                Y(order,i,j)=yorderi-y(order,j);

                ans=0.;

                //A-Jet calculation:
                if(parity==0){

                    for (int k=0; k<=evenOrder; k++) {
                        ans=ans+X(order-k,i,j)*X(k,i,j)+Y(order-k,i,j)*Y(k,i,j);
                    };

                    ans=2.*ans;
                    ans=ans+pow(X(order/2,i,j),2.)+pow(Y(order/2,i,j),2.);

                }
                else if(parity==1 || parity==-1){

                    for (int k=0; k<=oddOrder; k++) {
                        ans=ans+X(order-k,i,j)*X(k,i,j)+Y(order-k,i,j)*Y(k,i,j);
                    };

                    ans=2.*ans;

                }
                else std::cout << "ERROR: fillUpToOrder ERROR 1" << std::endl;

                A(order,i,j)=ans;

                //(B;DX,DY when i==0);C -Jets calculation:
                if(i==index_sa || j==index_sa){

                    dns=A(0,i,j);

                    if (order==0){

                        B(order,i,j)=pow(dns,-1.5);
                        C(order,i,j)=pow(dns,-2.5);

                    }
                    else if (order>0){

                        ans=0.;
                        bns=0.;

                        for (int k=0; k<=order-1; k++){

                            ans=ans+(-1.5*order+0.5*k)*A(order-k,i,j)*B(k,i,j);
                            bns=bns+(-2.5*order+1.5*k)*A(order-k,i,j)*C(k,i,j);

                        };

                        cns=order*dns;
                        B(order,i,j)=ans/cns;
                        C(order,i,j)=bns/cns;

                    }
                    else std::cout << "ERROR: fillUpToOrder ERROR 2" << std::endl;

                    ans=0.;
                    bns=0.;
                    for (int k=0; k<=order; k++){

                        cns=B(order-k,i,j)+Lambda*C(order-k,i,j);

                        ans=ans+X(k,i,j)*cns;
                        bns=bns+Y(k,i,j)*cns;

                    };

                    DX(order,i,j)= ans;
                    DY(order,i,j)= bns;

                }//end, if i==0 true
                else if(i!=index_sa || j!=index_sa){ //(B;DX,DY when i!=0)-Jets calculation

                    dns=A(0,i,j);

                    if (order==0) B(order,i,j)=pow(dns,-1.5);
                    else if (order>0){
                        ans=0.;
                        for (int k=0; k<=order-1; k++) ans=ans+(-1.5*order+0.5*k)*A(order-k,i,j)*B(k,i,j);
                        B(order,i,j)=ans/(order*dns);
                    }
                    else std::cout << "ERROR: fillUpToOrder ERROR 3" << std::endl;

                    ans=0.;
                    bns=0.;
                    for (int k=0; k<=order; k++){

                        cns=B(order-k,i,j);

                        ans=ans+X(k,i,j)*cns;
                        bns=bns+Y(k,i,j)*cns;

                    };

                    DX(order,i,j)= ans;
                    DY(order,i,j)= bns;

                };//end, if i!=0 true

            };//end, for j (j>i by construction)

        };//end, for i

        double factor=G/((order+2.)*(order+1.));
        int orderplustwo=order+2;

        x(orderplustwo,4) = ( m(0)*DX(order,0,4)+m(1)*DX(order,1,4)+m(2)*DX(order,2,4)+m(3)*DX(order,3,4) )*factor;
        y(orderplustwo,4) = ( m(0)*DY(order,0,4)+m(1)*DY(order,1,4)+m(2)*DY(order,2,4)+m(3)*DY(order,3,4) )*factor;

        x(orderplustwo,3) = ( m(0)*DX(order,0,3)+m(1)*DX(order,1,3)+m(2)*DX(order,2,3)                    )*factor;
        y(orderplustwo,3) = ( m(0)*DY(order,0,3)+m(1)*DY(order,1,3)+m(2)*DY(order,2,3)                    )*factor;

        x(orderplustwo,2) = ( m(0)*DX(order,0,2)+m(1)*DX(order,1,2)                   -m(3)*DX(order,2,3) )*factor;
        y(orderplustwo,2) = ( m(0)*DY(order,0,2)+m(1)*DY(order,1,2)                   -m(3)*DY(order,2,3) )*factor;

        x(orderplustwo,1) = ( m(0)*DX(order,0,1)                   -m(2)*DX(order,1,2)-m(3)*DX(order,1,3) )*factor;
        y(orderplustwo,1) = ( m(0)*DY(order,0,1)                   -m(2)*DY(order,1,2)-m(3)*DY(order,1,3) )*factor;

        x(orderplustwo,0) = ( -m(1)*x(orderplustwo,1)-m(2)*x(orderplustwo,2)-m(3)*x(orderplustwo,3) )/m(0);
        y(orderplustwo,0) = ( -m(1)*y(orderplustwo,1)-m(2)*y(orderplustwo,2)-m(3)*y(orderplustwo,3) )/m(0);

//        x(orderplustwo,0)=0.;
//        y(orderplustwo,0)=0.;

    };//end, for order

    TimeStepControl(n,N,x,y,thisIsTrue,infNorm,absoluteError,relativeError,timeStepControlMethod,delta_t);

    if( (t+delta_t)/T_Pr>=maxperiods ){
        //std::cout << "ALMOST DONE!" << " t/T=" << t/T_Pr << std::endl;
        double mymacheps=std::numeric_limits<double>::epsilon();
        delta_t=((maxperiods*T_Pr)-t)*(1.0000001);

    }

    TaylorSeriesExpansion(n,N,x,y,m,delta_t,xnew,ynew);

    t=t+delta_t;

}//end, TaylorMethod
