#include <iostream>
#include <cmath>
#include "matrix.h"
#include "cholesky.h"
#include "Sampler.h"
#include <ctime>

using namespace std;

int main() {

    srand(time(0));


    /////////////////////// PROBLEM 5.2 ///////////////////////
    const double sigma1=0.1,sigma2=0.2,sigma3=0.3,rho12=0.2,rho13=0.8,rho23=0.5;
    matrix A(3,3);
    A(0,0)=pow(sigma1,2);A(0,1)=sigma1*sigma2*rho12;A(0,2)=sigma1*sigma3*rho13;
    A(1,0)=A(0,1);A(1,1)=pow(sigma2,2);A(1,2)=sigma2*sigma3*rho23;
    A(2,0)=A(0,2);A(2,1)=A(1,2);A(2,2)=pow(sigma3,2);
    matrix L = cholesky(A);
    cout<<"Problem 5.2: Computing the Cholesky factorization"<<'\n'<<'\n';
    cout<<"L="<<'\n';
    print(L);
    cout<<'\n'<<'\n'<<"End of Problem 5.2"<<'\n'<<'\n';
    /////////////////////// PROBLEM 5.2 ///////////////////////



    /////////////////////// PROBLEM 6.3 ///////////////////////
    const int d=3;
    double k=10, T=1, r=0.01, K=100;
    double tk=T/k;

    NormalSampler stdnormal(0.0, 1.0);

    matrix sigma(d,1), S(k,d), a(d,1);
    sigma(-1,1)=sigma1; sigma(0,1)=sigma2; sigma(1,1)=sigma3;
    S(0,0)=100; S(0,1)=110; S(0,2)=120;// stock prices initial values
    a(-1,1)=0.333333333,a(0,1)=0.333333333,a(1,1)=0.333333333;

    for (int i=1;i<k;i++){
        for (int j=0;j<d;j++){
            double sum=0;
            for (int l=0;l<d;l++){
                sum=sum+L(j,l)*stdnormal.getnumber();
            }
            S(i,j)= S(i-1,j)*exp(  (r-0.5*pow(sigma(j-1,1),2))*tk + sum*sqrt(tk) );
        }
    }

    cout<<"Problem 6.3: Generating sample paths for a d-dimensional stochastic process "<<'\n'<<'\n';
    cout<<"S="<<'\n';
    print(S);
    cout<<'\n'<<'\n'<<"End of Problem 6.3"<<'\n'<<'\n';
    /////////////////////// PROBLEM 6.3 ///////////////////////

    /////////////////////// PROBLEM 7.3 ///////////////////////
    double V=0;
    double n1=100000;
    matrix v(n1,1);

    for (int kk=0;kk<n1;kk++){
        double fun=0;
        for (int i=0;i<d;i++){
            double sum1=0;
            for (int l=0;l<d;l++){
                sum1=sum1+L(i,l)*stdnormal.getnumber();
            }
            fun=fun+(a(i-1,1)*(S(0,i)*exp(  (r-(0.5*pow(sigma(i-1,1),2)))*T + sum1*sqrt(T) ) ) );
        }

        v(kk-1,1)=exp(-r*T)*fmax(fun-K,0);
        V=V+exp(-r*T)*fmax(fun-K,0);

    }
    V=V/n1;

    double sum2=0,s;
    for (int kk=0;kk<n1;kk++){
        sum2=sum2+pow(v(kk-1,1)-V,2);
    }
    s=sqrt( (1/(n1-1))*sum2  );

    double Vmin=V-((1.96*s)/sqrt(n1)), Vmax=V+((1.96*s)/sqrt(n1));

    cout<<"Problem 7.3: Computing the Basket call option price with a standard MC method "<<'\n'<<'\n';
    cout<<"Price V="<<V<<'\n';
    cout<<"95% Confidence interval:"<<"   "<<"("<<Vmin<<","<<Vmax<<")"<<'\n';
    cout<<s/sqrt(n1)<<'\n';
    cout<<'\n'<<'\n'<<"End of Problem 7.3"<<'\n'<<'\n';

    /////////////////////// PROBLEM 7.3 ///////////////////////

     /////////////////////// PROBLEM 7.4 ///////////////////////
    double Vc=0,Vc1=0,Vc2=0;
    double n2=100000;
    matrix vc(n2,1),vc1(n2,1),vc2(n2,1);
    double fun2=0,fun3=0,fun4=0;
    for (int kk=0;kk<n2;kk++){
        double funn=0;
        for (int i=0;i<d;i++){
            double sum1=0;
            for (int l=0;l<d;l++){
                sum1=sum1+L(i,l)*stdnormal.getnumber();
            }
            funn=funn+(a(i-1,1)*(S(0,i)*exp(  (r-(0.5*pow(sigma(i-1,1),2)))*T + sum1*sqrt(T) ) ) );
            if (i==0){
                fun2=(a(i-1,1)*(S(0,i)*exp(  (r-(0.5*pow(sigma(i-1,1),2)))*T + sum1*sqrt(T) ) ) )     +(a(0,1)*S(0,1))+(a(1,1)*S(0,2))-K;
            }

            else if (i==1){
                fun3=(a(i-1,1)*(S(0,i)*exp(  (r-(0.5*pow(sigma(i-1,1),2)))*T + sum1*sqrt(T) ) ) )     +(a(-1,1)*S(0,0))+(a(1,1)*S(0,2))-K;
            }
            else {
                fun4=(a(i-1,1)*(S(0,i)*exp(  (r-(0.5*pow(sigma(i-1,1),2)))*T + sum1*sqrt(T) ) ) )     +(a(-1,1)*S(0,0))+(a(0,1)*S(0,1))-K;
            }

        }
        double beta=1;
        vc(kk-1,1)=(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun2,0) - 10.224);
        Vc=Vc+(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun2,0) - 10.224);

        vc1(kk-1,1)=(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun3,0) - 10.389);
        Vc1=Vc1+(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun3,0) - 10.389);

        vc2(kk-1,1)=(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun4,0) - 11.161);
        Vc2=Vc2+(exp(-r*T)*fmax(funn-K,0)) - beta*(exp(-r*T)*fmax(fun4,0) - 11.161);
    }
    Vc=Vc/n2;
    Vc1=Vc1/n2;
    Vc2=Vc2/n2;

    double sum3=0,sum4=0,sum5=0,ss,ss1,ss2;
    for (int kk=0;kk<n2;kk++){
        sum3=sum3+pow(vc(kk-1,1)-Vc,2);
        sum4=sum4+pow(vc1(kk-1,1)-Vc1,2);
        sum5=sum5+pow(vc2(kk-1,1)-Vc2,2);
    }
    ss=sqrt( (1/(n2-1))*sum3  );
    ss1=sqrt( (1/(n2-1))*sum4  );
    ss2=sqrt( (1/(n2-1))*sum5  );

    double Vcmin=Vc-((1.96*ss)/sqrt(n2)), Vcmax=Vc+((1.96*ss)/sqrt(n2));
    double Vc1min=Vc1-((1.96*ss1)/sqrt(n2)), Vc1max=Vc1+((1.96*ss1)/sqrt(n2));
    double Vc2min=Vc2-((1.96*ss2)/sqrt(n2)), Vc2max=Vc2+((1.96*ss2)/sqrt(n2));

    cout<<"Problem 7.4: Computing the Basket call option price using control variates "<<'\n'<<'\n';
    cout<<"Price Vc1="<<Vc<<'\n';
    cout<<"95% Confidence interval:"<<"   "<<"("<<Vcmin<<","<<Vcmax<<")"<<'\n';
    cout<<"Standard Error:"<<"   "<<ss/sqrt(n2)<<'\n';
    cout<<"Price Vc2="<<Vc1<<'\n';
    cout<<"95% Confidence interval:"<<"   "<<"("<<Vc1min<<","<<Vc1max<<")"<<'\n';
    cout<<"Standard Error:"<<"   "<<ss1/sqrt(n2)<<'\n';
    cout<<"Price Vc3="<<Vc2<<'\n';
    cout<<"95% Confidence interval:"<<"   "<<"("<<Vc2min<<","<<Vc2max<<")"<<'\n';
    cout<<"Standard Error:"<<"   "<<ss2/sqrt(n2)<<'\n';
    cout<<'\n'<<'\n'<<"End of Problem 7.4"<<'\n'<<'\n';

    ///////////////////// PROBLEM 7.4 ///////////////////////


    return 0;
}
