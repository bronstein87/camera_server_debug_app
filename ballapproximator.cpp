#include "ballapproximator.h"
#include <math.h>
BallApproximator::BallApproximator()
{

}

void BallApproximator::readInitData(const QString &path, double deltaT,bool writeToFile)
{
    double lenght;
    int I,J;
    FILE *FL;

    FL=fopen(path.toStdString().c_str(), "r");
    fscanf(FL,"%lf %lf %lf", &CAM0[0], &CAM0[1], &CAM0[2]);
    fscanf(FL,"%lf %lf %lf", &CAM1[0], &CAM1[1], &CAM1[2]);
    if(writeToFile){
        fprintf(FLIST,"coordiates of cameras: c0 & c1\n");
    }
    //fprintf(FLIST,"coordiates of cameras: c0 & c1\n");
    for (I=0; I<3; I++) //fprintf(FLIST,"%d %9.5f %9.5f\n",I,CAM0[I],CAM1[I]);

        if(writeToFile){
            fprintf(FLIST,"%d %9.5f %9.5f\n",I,CAM0[I],CAM1[I]);
            fprintf(FLIST,"\ndata from c0\n");
        }
    //fprintf(FLIST,"\ndata from c0\n");
    fscanf(FL,"%d",&N0);

    if(writeToFile){
        fprintf(FLIST,"%d\n",N0);
    }
    // fprintf(FLIST,"%d\n",N0);
    for (I=0; I<N0; I++) fscanf(FL,"%lf %lf %lf",&U0[I][0],&U0[I][1],&U0[I][2]);
    for (I=0; I<N0; I++)
        if(writeToFile){
            fprintf(FLIST,"%2d %11.7f %11.7f %11.7f\n",I,U0[I][0],U0[I][1],U0[I][2]);
            fprintf(FLIST,"\ndata from c1\n");
            fprintf(FLIST,"%d\n",N1);
        }
    // fprintf(FLIST,"%2d %11.7f %11.7f %11.7f\n",I,U0[I][0],U0[I][1],U0[I][2]);

    //fprintf(FLIST,"\ndata from c1\n");
    fscanf(FL,"%d",&N1);
    //fprintf(FLIST,"%d\n",N1);
    for (I=0; I<N1; I++) fscanf(FL,"%lf %lf %lf",&U1[I][0],&U1[I][1],&U1[I][2]);
    for (I=0; I<N1; I++)
        if(writeToFile){
            fprintf(FLIST,"%2d %11.7f %11.7f %11.7f\n",I,U1[I][0],U1[I][1],U1[I][2]);
        }
    for (I=0; I<N0; I++)
    {
        T0[I]=I * deltaT;
        lenght = 0;

        for (J=0; J<3; J++)
        {
            lenght += U0[I][J] * U0[I][J];
        }

        lenght = sqrt(lenght);

        for (J=0; J<3; J++)
        {
            U0[I][J] /= lenght;
        }
    }

    for (I=0; I<N1; I++)
    {
        T1[I]=I * deltaT;
        lenght = 0;
        for (J=0; J<3; J++) lenght+=U1[I][J]*U1[I][J];
        lenght=sqrt(lenght);
        for (J=0; J<3; J++) U1[I][J]/=lenght;
    }
}
#include <QDebug>
void BallApproximator::readData(QVector<Calibration::Position> f, QVector<Calibration::Position> s,
                                QVector <double> fTime, QVector <double> sTime,
                                Calibration::Position fPos, Calibration::Position sPos)
{

    fR = f;
    sR = s;
    fTimeR = fTime;
    sTimeR = sTime;
    fPosR = fPos;
    sPosR = sPos;


    N0 = f.size();
    N1 = s.size();
    CAM0[0] = fPos.X;
    CAM0[1] = fPos.Y;
    CAM0[2] = fPos.Z;

    CAM1[0] = sPos.X;
    CAM1[1] = sPos.Y;
    CAM1[2] = sPos.Z;
    for (qint32 i = 0; i < f.size(); ++i)
    {
        U0[i][0] = f[i].X;
        U0[i][1] = f[i].Y;
        U0[i][2] = f[i].Z;
        T0[i] = fTime[i];
    }

    for (qint32 i = 0; i < s.size(); ++i)
    {
        U1[i][0] = s[i].X;
        U1[i][1] = s[i].Y;
        U1[i][2] = s[i].Z;
        T1[i] = sTime[i];
    }


    for (qint32 i = 0; i < f.size(); i++)
    {

        double lenght = 0;
        for (qint32 j = 0; j < 3; j++)
        {
            lenght += U0[i][j] * U0[i][j];
        }

        lenght = sqrt(lenght);

        for (qint32 j = 0; j < 3; j++)
        {
            U0[i][j] /= lenght;
        }
    }

    for (qint32 i = 0; i < s.size(); i++)
    {
        double lenght = 0;
        for (qint32 j = 0; j < 3; j++)
        {
            lenght += U1[i][j] * U1[i][j];
        }

        lenght = sqrt(lenght);

        for (qint32 j = 0; j < 3; j++)
        {
            U1[i][j] /= lenght;
        }
    }
    qDebug() << CAM0[0] << CAM0[1] << CAM0[2];
    for (int i = 0; i < f.size(); ++i)
    {
        qDebug() << U0[i][0] << U0[i][1] << U0[i][2] << QString::number(T0[i],'g', 11);
    }
    qDebug() << "///////////";
    qDebug() << CAM1[0] << CAM1[1] << CAM1[2];
    for (int j = 0; j < s.size(); j++)
    {
        qDebug() << U1[j][0] << U1[j][1] << U1[j][2] << QString::number(T1[j],'g', 11);
    }

    TIN=T0[N0-1];
    if (T1[N1-1]<TIN) TIN=T1[N1-1];
    fprintf(FLIST,"TIN= %11.7f\n",TIN);
}

void BallApproximator::rotateMovementParameters(double pos[], double v[], double a[])
{
    double posTmp[3], vTmp[3], aTmp[3];
    posTmp[0] = xNonLinearParams[0];
    vTmp[0] = xNonLinearParams[1];
    aTmp[0] = xNonLinearParams[2];

    posTmp[1] = yNonLinearParams[0];
    vTmp[1] = yNonLinearParams[1];
    aTmp[1] = yNonLinearParams[2];

    posTmp[2] = zNonLinearParams[0];
    vTmp[2] = zNonLinearParams[1];
    aTmp[2] = zNonLinearParams[2];
    double mInit [3][3] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double mRot[3][3];
    rotateOZ(45.0 * degreesToRad, mInit, mRot);

    multMatrixVector(mRot, posTmp, pos);
    multMatrixVector(mRot, vTmp, v);
    multMatrixVector(mRot, aTmp, a);
}

//void BallApproximator::calculateIntersect(double distFirst, double distSecond, double &t1, double &t2)
//{

//}




void BallApproximator::DIST(double *U, double *V, double *POS, double *E)
{
    int I;
    double UV,X,Y,P,Q,P0[3],P1[3];

    UV=0; X=0; Y=0;
    for (I=0; I<3; I++) {
        UV+=U[I]*V[I]; X+=DC[I]*U[I]; Y+=DC[I]*V[I];
    }
    P=(X-UV*Y)/(1-UV*UV); Q=-(Y-UV*X)/(1-UV*UV);
    Y=0;
    for (I=0; I<3; I++) {
        P0[I]=CAM0[I]+P*U[I]; P1[I]=CAM1[I]+Q*V[I];
        POS[I]=0.5*(P0[I]+P1[I]); X=P0[I]-P1[I]; Y+=X*X;
    }
    *E=sqrt(Y);
}



//void BallApproximator::PARAB(int K, double *T, double *X, double *RVW, double *E)
//{
//    int I;
//    double S0,S1,S2,S3,S4,Y0,Y1,Y2,A,B,D,A0,A1,A2,L0,L1,L2,L3,L4,L5;

//    S0=K; S1=0; S2=0; S3=0; S4=0; A0=0; A1=0; A2=0;
//    for (I=0; I<K; I++) {
//        A=T[I] - TIN; B=X[I]; D=A*A; S1+=A; S2+=D; S3+=A*D; S4+=D*D;
//        A0+=B; A1+=A*B; A2+=B*D;
//    }
//    L0=sqrt(S0); L1=S1/L0; L2=S2/L0;
//    L3=sqrt(S2-L1*L1); L4=(S3-L1*L2)/L3; L5=sqrt(S4-L2*L2-L4*L4);
//    Y0=A0/L0; Y1=(A1-L1*Y0)/L3; Y2=(A2-L2*Y0-L4*Y1)/L5;
//    A2=Y2/L5; A1=(Y1-L4*A2)/L3; A0=(Y0-L1*A1-L2*A2)/L0;
//    B=0;
//    for (I=0; I<K; I++) {A=X[I]-A0-T[I]*(A1+T[I]*A2); B+=A*A;}
//    *E=sqrt(B/(K-3));
//    RVW[0]=A0; RVW[1]=A1; RVW[2]=A2;
//}

void BallApproximator::PARAB(int K, double *T, double *X, double *RVW, double *E)
{
    int I;
    double S0,S1,S2,S3,S4,Y0,Y1,Y2,A,B,D,A0,A1,A2,L0,L1,L2,L3,L4,L5;

    for (I=0; I<K; I++) fprintf(FLIST,"%2d  %12.6f  %12.6f\n",I,T[I],X[I]);
    //getchar();
    S0=K; S1=0; S2=0; S3=0; S4=0; A0=0; A1=0; A2=0;
    for (I=0; I<K; I++) {
        A=T[I]-TIN; B=X[I]; D=A*A; S1+=A; S2+=D; S3+=A*D; S4+=D*D;
        A0+=B; A1+=A*B; A2+=B*D;
    }
    L0=sqrt(S0); L1=S1/L0; L2=S2/L0;
    L3=sqrt(S2-L1*L1); L4=(S3-L1*L2)/L3; L5=sqrt(S4-L2*L2-L4*L4);
    Y0=A0/L0; Y1=(A1-L1*Y0)/L3; Y2=(A2-L2*Y0-L4*Y1)/L5;
    A2=Y2/L5; A1=(Y1-L4*A2)/L3; A0=(Y0-L1*A1-L2*A2)/L0;
    B=0;
    for (I=0; I<K; I++) {
        D=T[I]-TIN; A=X[I]-A0-D*(A1+D*A2); B+=A*A;
    }
    *E=sqrt(B/(K-3)); RVW[0]=A0; RVW[1]=A1; RVW[2]=A2;
    fprintf(FLIST,"x,vx,ax/2,rms: %12.5e %12.5e %12.5e %10.3e\n",A0,A1,A2,*E);
    for (I=0; I<K; I++) {
        D=T[I]-TIN; S0=A0+D*(A1+D*A2); S1=A1+2*A2*D;
        fprintf(FLIST,"%2d  %12.5f  %12.5e  %12.5e  %12.5e\n",I,D,S0,S1,X[I]-S0);
    }
    //getchar();
}



int BallApproximator::CHOLDET1(int N, double (*A)[nonLinearParamsCount], double *P)
{
    int I,J,K;
    double X;

    for (I=0; I<N; I++)
        for (J=I; J<N; J++) {
            X=A[I][J];
            for (K=0; K<I; K++) X=X-A[J][K]*A[I][K];
            if (J==I) {
                if (X<=0.0) return 0;
                P[I]=1/sqrt(X);
            }
            else A[J][I]=X*P[I];
        }
    return 1;
}

void BallApproximator::CHOLSOL1(int N, int R, double (*A)[nonLinearParamsCount], double *P, double (*B)[numberResults], double (*X)[numberResults])
{
    int I,J,K;
    double Z;

    for (J=0; J<R; J++) {
        //solution of LY==B
        for (I=0; I<N; I++) {
            Z=B[I][J];
            for (K=0; K<I; K++) Z=Z-A[I][K]*X[K][J];
            X[I][J]=Z*P[I];
        }
        //solution of UX==Y
        for (I=N-1; I>=0; I--) {
            Z=X[I][J];
            for (K=I+1; K<N; K++) Z=Z-A[K][I]*X[K][J];
            X[I][J]=Z*P[I];
        }
    }
}

int BallApproximator::CHOLINV1(int N, double (*A)[nonLinearParamsCount])
{
    int I,J,K,I1,J1;
    double X,Y,Z;

    //formation of L
    for (I=0; I<N; I++) {
        I1=I+1;
        for (J=I; J<N; J++) {
            J1=J+1; X=A[I][J];
            for (K=I-1; K>=0; K--) X=X-A[J1][K]*A[I1][K];
            if (J==I) {
                if (X<=0) return 0;
                Y=1/sqrt(X); A[I1][I]=Y;
            }
            else A[J1][I] = X * Y;
        }
    }

    //inversion of L
    for (I=0; I<N; I++) for (J=I+1; J<N; J++) {
        Z=0; J1=J+1;
        for (K=J-1; K>=I; K--) Z=Z-A[J1][K]*A[K+1][I];
        A[J1][I]=Z*A[J1][J];
    }

    //calculation of the inverse of A
    for (I=0; I<N; I++) for (J=0; J<N; J++) {
        Z=0; J1=N+1;
        for (K=J+1; K<J1; K++) Z=Z+A[K][J]*A[K][I];
        A[J+1][I]=Z;
    }
    return 1;
}

int BallApproximator::JACOBI(int N, int EIVEC, double (*A)[nonLinearParamsCount], double (*V)[nonLinearParamsCount], double *D)
{
    double SM,C,S,T,H,G,TAU,THETA,TRESH;
    int I,J,P,Q,ROT;
    double B[nonLinearParamsCount],Z[nonLinearParamsCount];

    if (EIVEC==1) {
        for (P=0; P<N; P++)
            for (Q=0; Q<N; Q++) {
                if (P==Q) V[P][Q]=1; else V[P][Q]=0;
            }
    }
    for (P=0; P<N; P++) {
        B[P]=A[P][P]; D[P]=B[P]; Z[P]=0;
    }
    ROT=0;
    for (I=0; I<50; I++) { //swp
        SM=0;
        for (P=0; P<N-1; P++)
            for (Q=P+1; Q<N; Q++) SM=SM+fabs(A[P][Q]);
        if (qFuzzyCompare(SM, 0.0)) return ROT;
        if (I<4) TRESH=0.2*SM/N/N; else TRESH=0;
        for (P=0; P<N-1; P++)
            for (Q=P+1; Q<N; Q++) {
                G=100*fabs(A[P][Q]);
                if (I>4 && fabs(D[P])+G==fabs(D[P]) && fabs(D[Q])+G==fabs(D[Q]))
                    A[P][Q]=0;
                else
                    if (fabs(A[P][Q])>TRESH) { //rotate
                        H=D[Q]-D[P];
                        if (fabs(H)+G==fabs(H)) T=A[P][Q]/H;
                        else {
                            THETA=0.5*H/A[P][Q];
                            T=1/(fabs(THETA)+sqrt(1+THETA*THETA));
                            if (THETA<0) T=-T;
                        } //calculation of tan of rotataion angle
                        C=1/sqrt(1+T*T); S=T*C;
                        TAU=S/(1+C); H=T*A[P][Q];
                        Z[P]=Z[P]-H; Z[Q]=Z[Q]+H;
                        D[P]=D[P]-H; D[Q]=D[Q]+H;
                        A[P][Q]=0;
                        for (J=0; J<=P-1; J++) {
                            G=A[J][P]; H=A[J][Q];
                            A[J][P]=G-S*(H+G*TAU);
                            A[J][Q]=H+S*(G-H*TAU);
                        } //1<=J<P
                        for (J=P+1; J<=Q-1; J++) {
                            G=A[P][J]; H=A[J][Q];
                            A[P][J]=G-S*(H+G*TAU);
                            A[J][Q]=H+S*(G-H*TAU);
                        } //P<J<Q
                        for (J=Q+1; J<N; J++) {
                            G=A[P][J]; H=A[Q][J];
                            A[P][J]=G-S*(H+G*TAU);
                            A[Q][J]=H+S*(G-H*TAU);
                        } //Q<J<=N
                        if (EIVEC==1)
                            for (J=0; J<N; J++) {
                                G=V[J][P];  H=V[J][Q];
                                V[J][P]=G-S*(H+G*TAU);
                                V[J][Q]=H+S*(G-H*TAU);
                            } //calculation of V
                        ROT=ROT+1;
                    } //rotate
            }
        for (P=0; P<N; P++) {
            B[P]=B[P]+Z[P]; D[P]=B[P]; Z[P]=0;
        } //p
    } //swp
    return ROT;
}


int BallApproximator::CHOLDET2(int N, double *A)
{
    int I,J,K,P,Q,R;
    double X;

    P=-1;
    for (I=0; I<N; I++) {
        Q=P+1; R=-1;
        for (J=0; J<=I; J++) {
            X=A[P+1];
            for (K=Q; K<=P; K++) {R++; X-=A[K]*A[R];}
            R++; P++;
            if (J==I) {
                if (X<=0) return 0;
                A[P]=1/sqrt(X);
            }
            else A[P]=X*A[R];
        }
    }
    return 1;
}

void BallApproximator::CHOLSOL2(int N, int R, double *A, double (*B)[numberResults])
{
    int I,J,K,P,S;
    double X;

    for (J=0; J<R; J++) {
        //solution of LY=B
        P=0;
        for (I=0; I<N; I++) {
            X=B[I][J];
            for (K=0; K<I; K++) {X-=A[P]*B[K][J]; P++;}
            B[I][J]=X*A[P]; P++;
        }

        //solution of UX=Y
        for (I=N-1; I>=0; I--) {
            P--; S=P; X=B[I][J];
            for (K=N-1; K>I; K--) {X-=A[S]*B[K][J]; S=S-K;}
            B[I][J]=X*A[S];
        }
    }
}

int BallApproximator::CHOLINV2(int N, double *A)
{
    int I,J,K,P,Q,R,S,T,U;
    double X,Y;

    //formation of L
    P=-1;
    for (I=0; I<N; I++) {
        Q=P+1; R=-1;
        for (J=0; J<=I; J++) {
            X=A[P+1];
            for (K=Q; K<=P; K++) {R++; X-=A[K]*A[R];}
            R++; P++;
            if (J==I) {
                if (X<=0) return 0;
                A[P]=1/sqrt(X);
            }
            else A[P]=X*A[R];
        }
    }

    //inversion of L
    P=-1; R=-1; T=-1;
    for (I=1; I<N; I++) {
        P++; R=R+I+1;
        Y=A[R+1];
        for (J=1; J<=I; J++) {
            P++; T++; S=T;
            U=I-1; X=0;
            for (K=R; K>=P; K--) {
                X-=A[K]*A[S];
                S-=U; U--;
            }
            A[P]=X*Y;
        }
    }

    //calculation of the inverse of A
    P=-1;
    for (I=0; I<N; I++) {
        R=P+I+1;
        for (J=0; J<=I; J++) {
            S=R; P++; T=P;
            X=0;
            for (K=I; K<N; K++) {
                X+=A[S]*A[T];
                S+=K+1; T+=K+1;
            }
            A[P]=X;
        }
    }

    return 1;
}

//void BallApproximator::NORMEQ(double *X, double (*Z)[nonLinearParamsCount], double (*V)[1], double *F)
//{
//    double A[3][nonLinearParamsCount],P[3] ,S,RB, DT;
//    int I,J,K,L;

//    *F=0;
//    for (I=0; I<nonLinearParamsCount; I++) {
//        V[I][0]=0; for (J=I; J<nonLinearParamsCount; J++) Z[I][J]=0;
//    }

//    for (I=0; I<N0; I++) {
//        S=0;DT=T0[I]-TIN;
//        for (J=0; J<3; J++) {
//            P[J]=X[J]+T0[I]*(X[J+3]+T0[I]*X[J+6])-CAM0[J];
//            S+=P[J]*P[J];
//        }
//        S=sqrt(S);
//        for (J=0; J<3; J++) {P[J]/=S; RB=U0[I][J]-P[J]; *F+=RB*RB;}
//        for (J=0; J<3; J++) for (K=J; K<3; K++) {
//            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
//        }
//        for (J=0; J<3; J++) for (K=0; K<3; K++) {
//            A[J][K+3]=T0[I]*A[J][K]; A[J][K+6]=T0[I]*A[J][K+3];
//        }
//        for (J=0; J<nonLinearParamsCount; J++) {
//            for (L=0; L<3; L++) V[J][0]+=A[L][J]*(U0[I][L]-P[L]);
//            for (K=J; K<nonLinearParamsCount; K++) for (L=0; L<3; L++) Z[J][K]+=A[L][J]*A[L][K];
//        }
//    }

//    for (I=0; I<N1; I++) {
//        S=0;DT=T1[I]-TIN;
//        for (J=0; J<3; J++) {
//            P[J]=X[J]+T1[I]*(X[J+3]+T1[I]*X[J+6])-CAM1[J];
//            S+=P[J]*P[J];
//        }
//        S=sqrt(S);
//        for (J=0; J<3; J++) {P[J]/=S; RB=U1[I][J]-P[J]; *F+=RB*RB;}
//        for (J=0; J<3; J++) for (K=J; K<3; K++) {
//            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
//        }
//        for (J=0; J<3; J++) for (K=0; K<3; K++) {
//            A[J][K+3]=T1[I]*A[J][K]; A[J][K+6]=T1[I]*A[J][K+3];
//        }
//        for (J=0; J<nonLinearParamsCount; J++) {
//            for (L=0; L<3; L++) V[J][0]+=A[L][J]*(U1[I][L]-P[L]);
//            for (K=J; K<nonLinearParamsCount; K++) for (L=0; L<3; L++) Z[J][K]+=A[L][J]*A[L][K];
//        }
//    }

//}




void BallApproximator::NORMEQ(double *X, double (*Z)[nonLinearParamsCount], double (*V)[1], double *F)
{
    double A[3][nonLinearParamsCount],P[3],S,RB,DT;
    int I,J,K,L;

    *F=0;
    for (I=0; I<nonLinearParamsCount; I++) {
        V[I][0]=0; for (J=I; J<nonLinearParamsCount; J++) Z[I][J]=0;
    }

    for (I=0; I<N0; I++) {
        S=0; DT=T0[I]-TIN;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*(X[J+3]+DT*X[J+6])-CAM0[J]; S+=P[J]*P[J];
        }
        S=sqrt(S);
        for (J=0; J<3; J++) {P[J]/=S; RB=U0[I][J]-P[J]; *F+=RB*RB;}
        for (J=0; J<3; J++) for (K=J; K<3; K++) {
            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
        }
        for (J=0; J<3; J++) for (K=0; K<3; K++) {
            A[J][K+3]=DT*A[J][K]; A[J][K+6]=DT*A[J][K+3];
        }
        for (J=0; J<nonLinearParamsCount; J++) {
            for (L=0; L<3; L++) V[J][0]+=A[L][J]*(U0[I][L]-P[L]);
            for (K=J; K<nonLinearParamsCount; K++) for (L=0; L<3; L++) Z[J][K]+=A[L][J]*A[L][K];
        }
    }

    for (I=0; I<N1; I++) {
        S=0; DT=T1[I]-TIN;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*(X[J+3]+DT*X[J+6])-CAM1[J]; S+=P[J]*P[J];
        }
        S=sqrt(S);
        for (J=0; J<3; J++) {P[J]/=S; RB=U1[I][J]-P[J]; *F+=RB*RB;}
        for (J=0; J<3; J++) for (K=J; K<3; K++) {
            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
        }
        for (J=0; J<3; J++) for (K=0; K<3; K++) {
            A[J][K+3]=DT*A[J][K]; A[J][K+6]=DT*A[J][K+3];
        }
        for (J=0; J<nonLinearParamsCount; J++) {
            for (L=0; L<3; L++) V[J][0]+=A[L][J]*(U1[I][L]-P[L]);
            for (K=J; K<nonLinearParamsCount; K++) for (L=0; L<3; L++) Z[J][K]+=A[L][J]*A[L][K];
        }
    }

}



void BallApproximator::SOLVER(double *X,bool writeToFile)
{
    double EPSI=1e-6,W,NEV,F,
            Z[nonLinearParamsCount][nonLinearParamsCount],
            V[nonLinearParamsCount][1],
            P[nonLinearParamsCount];
    int NN,I,RES;

    NN=0;
    do {
        NN=NN+1;
        NORMEQ(X,Z,V,&F);

        {
            int J;
            double U[nonLinearParamsCount][nonLinearParamsCount];

            for(I=0; I<nonLinearParamsCount; I++) for(J=I; J<nonLinearParamsCount; J++) U[I][J]=Z[I][J];
            RES=JACOBI(nonLinearParamsCount,0,U,U,P);
            if(writeToFile){
                printf("eigen values (ROT=%d)\n",RES);}
            for (I=0; I<nonLinearParamsCount; I++)
            {
                if(writeToFile){
                    printf("%12.5e ",P[I]); printf("\n");}
            }

            if(writeToFile){
                fprintf(FLIST,"eigen values (ROT=%d)\n",RES);}
            for (I=0; I<nonLinearParamsCount; I++)
            {
                if(writeToFile){
                    fprintf(FLIST,"%12.5e\n",P[I]); fprintf(FLIST,"\n");}
            }
        }
        NEV=0;
        for (I=0; I<nonLinearParamsCount; I++) NEV+=V[I][0]*V[I][0];
        NEV=sqrt(NEV);
        RES=CHOLDET1(nonLinearParamsCount,Z,P);
        if (RES) {
            CHOLSOL1(nonLinearParamsCount,1,Z,P,V,V);
            W=0;
            for (I=0; I<nonLinearParamsCount; I++) {
                X[I]+=V[I][0]; W+=V[I][0]*V[I][0];
            }
            if(writeToFile){
                printf("N=%2d step=%14.7e nev=%14.7e F=%14.7e\n",NN,W,NEV,F);
                fprintf(FLIST,"N=%2d step=%14.7e nev=%14.7e F=%14.7e\n",NN,W,NEV,F);}
            for (I=0; I<nonLinearParamsCount; I++)
                if(writeToFile){
                    printf("X[%d]=%17.10e;\n",I,X[I]);
                    printf("\n");}
            for (I=0; I<nonLinearParamsCount; I++) //fprintf(FLIST,"X[%d]=%17.10e;\n",I,X[I]);
                if(writeToFile){
                    fprintf(FLIST,"X[%d]=%17.10e;\n",I,X[I]);
                    fprintf(FLIST,"\n");}
        } else
        {

            if(writeToFile){
                printf("RES=0\n");
            }
            getchar();
            return;
        }
    } while (W>EPSI);

    xNonLinearParams[0] = X[0];
    xNonLinearParams[1] = X[3];
    xNonLinearParams[2] = X[6];
    yNonLinearParams[0] = X[1];
    yNonLinearParams[1] = X[4];
    yNonLinearParams[2] = X[7];
    zNonLinearParams[0] = X[2];
    zNonLinearParams[1] = X[5];
    zNonLinearParams[2] = X[8];

}




void BallApproximator::COVMAT(int NPAR, int NDIV, double *X,bool writeToFile)
{
    int I,J,K,ROT,RES;
    double F,DH,Z[nonLinearParamsCount][nonLinearParamsCount],Z1[nonLinearParamsCount][nonLinearParamsCount],ZV[nonLinearParamsCount][nonLinearParamsCount],ZIN[nonLinearParamsCount+1][nonLinearParamsCount];
    double AD[nonLinearParamsCount],V1[nonLinearParamsCount][1];
    if(writeToFile){
        printf("\n");
        fprintf(FLIST,"\n");
        printf("COVMAT: dim=%d, degree of freedom=%d\n",NPAR,NDIV);
        fprintf(FLIST,"COVMAT: dim=%d, degree of freedom=%d\n",NPAR,NDIV);}
    NORMEQ(X,Z1,V1,&F);
    DH=sqrt(F/NDIV);
    if(writeToFile){
        printf("FCRIT %14.7e, DH %11.7e\n",F,DH);
        fprintf(FLIST,"FCRIT %14.7e, sigma %11.7e\n",F,DH);
        printf("\n");
        fprintf(FLIST,"\n");}

    for (I=0; I<NPAR; I++) for (J=I; J<NPAR; J++) {
        Z[I][J]=Z1[I][J]; ZIN[I][J]=Z[I][J];
    }
    ROT=JACOBI(NPAR,1,Z,ZV,AD);
    if(writeToFile){
        printf("eigen values and the number of rotations\n");}
    for (K=0; K<NPAR/6; K++) {
        for (I=6*K; I<6*K+6; I++) printf("%12.5e ",AD[I]);
        if(writeToFile){
            printf("\n");}
    }
    for (I=6*(NPAR/6); I<NPAR; I++) printf("%12.5e ",AD[I]);

    if(writeToFile){
        printf(" %d\n",ROT);
        printf("\n");

        printf("eigen vectors\n");}
    for (K=0; K<NPAR/6; K++) {
        for (I=6*K; I<6*K+6; I++)

            if(writeToFile){
                printf("%12.5e ",AD[I]);
                printf("\n");
                printf("\n");}

        for (I=0; I<NPAR; I++) {
            for (J=6*K; J<6*K+6; J++) printf("%12.5e ",ZV[I][J]);
            if(writeToFile){
                printf("\n");}
        }
        if(writeToFile){
            printf("\n");}
    }
    if (NPAR%6>0) {
        for (I=6*(NPAR/6); I<NPAR; I++) printf("%12.5e ",AD[I]);

        if(writeToFile){
            printf("\n");
            printf("\n");}
        for (I=0; I<NPAR; I++) {
            for (J=6*(NPAR/6); J<NPAR; J++) printf("%12.5e ",ZV[I][J]);
            if(writeToFile){
                printf("\n");}
        }
        if(writeToFile){
            printf("\n");}
    }
    if(writeToFile){
        fprintf(FLIST,"eigen values and the number of rotations\n");}
    for (K=0; K<NPAR/6; K++) {
        for (I=6*K; I<6*K+6; I++)

            if(writeToFile){
                fprintf(FLIST,"%12.5e ",AD[I]);
                fprintf(FLIST,"\n");}
    }
    for (I=6*(NPAR/6); I<NPAR; I++)
        if(writeToFile){
            fprintf(FLIST,"%12.5e ",AD[I]);
            fprintf(FLIST," %d\n",ROT);
            fprintf(FLIST,"\n");
            fprintf(FLIST,"eigen vectors\n");}
    for (K=0; K<(NPAR/6); K++) {
        for (I=6*K; I<6*K+6; I++)
            if(writeToFile){
                fprintf(FLIST,"%12.5e ",AD[I]);
                fprintf(FLIST,"\n");
                fprintf(FLIST,"\n");}

        for (I=0; I<NPAR; I++) {
            for (J=6*K; J<6*K+6; J++)
                if(writeToFile){
                    fprintf(FLIST,"%12.5e ",ZV[I][J]);
                    fprintf(FLIST,"\n");}
        }
        if(writeToFile){
            fprintf(FLIST,"\n");}
    }
    if (NPAR%6>0) {
        for (I=6*(NPAR/6); I<NPAR; I++)
            if(writeToFile){
                fprintf(FLIST,"%12.5e ",AD[I]);
                fprintf(FLIST,"\n");
                fprintf(FLIST,"\n");}
        for (I=0; I<NPAR; I++) {
            for (J=6*(NPAR/6); J<NPAR; J++)
                if(writeToFile){
                    fprintf(FLIST,"%12.5e ",ZV[I][J]);
                    fprintf(FLIST,"\n");}
        }
        if(writeToFile){
            fprintf(FLIST,"\n");}
    }
    RES=CHOLINV1(NPAR,ZIN);
    if(writeToFile){
        printf("RES %d\n",RES);}
    for (I=0; I<NPAR; I++) AD[I]=DH*sqrt(ZIN[I+1][I]);
    if(writeToFile){
        printf("standart deviations\n");}
    for (K=0; K<(NPAR/6); K++) {
        for (I=6*K; I<6*K+6; I++)
            if(writeToFile){
                printf("%12.5e ",AD[I]);
                printf("\n");}
    }
    for (I=6*(NPAR/6); I<NPAR; I++)
        if(writeToFile){
            printf("%12.5e ",AD[I]);
            printf("\n");
            fprintf(FLIST,"standart deviations\n");
        }
    for (K=0; K<(NPAR/6); K++) {
        for (I=6*K; I<6*K+6; I++)
            if(writeToFile){
                fprintf(FLIST,"%12.5e ",AD[I]);
                fprintf(FLIST,"\n");}
    }
    for (I=6*(NPAR/6); I<NPAR; I++)
        if(writeToFile)
        {
            fprintf(FLIST,"%12.5e ",AD[I]);
            fprintf(FLIST,"\n");
        }
}


void BallApproximator::FIRST()
{
    double TB,TF,Y[maxNumberOfMeasures],P[3],E;
    int I,J;

    TB=T0[0]; if(TB<T1[0]) TB=T1[0];
    TF=T0[N0-1]; if(TF>T1[N1-1])TF=T1[N1-1];
    E=(TF-TB)/(fisrtMeasureCount-1);
    for (I=0; I<fisrtMeasureCount; I++) timeSame[I]=TB+E*I;

    for (J=0; J<3; J++) {
        for (I=0; I<N0; I++) Y[I]=U0[I][J];
        PARAB(N0,T0,Y,P,&E);
        for (I=0; I<fisrtMeasureCount; I++) U0[I][J]=P[0]+(timeSame[I]-TIN)*(P[1]+(timeSame[I]-TIN)*P[2]);
        for (I=0; I<N1; I++) Y[I]=U1[I][J];
        PARAB(N1,T1,Y,P,&E);
        for (I=0; I<fisrtMeasureCount; I++) U1[I][J]=P[0]+(timeSame[I]-TIN)*(P[1]+(timeSame[I]-TIN)*P[2]);
    }

    for (I=0; I<fisrtMeasureCount; I++) {
        E=0;
        for (J=0; J<3; J++) E+=U0[I][J]*U0[I][J];
        E=sqrt(E);
        for (J=0; J<3; J++) U0[I][J]/=E;
    }
    for (I=0; I<fisrtMeasureCount; I++) {
        E=0;
        for (J=0; J<3; J++) E+=U1[I][J]*U1[I][J];
        E=sqrt(E);
        for (J=0; J<3; J++) U1[I][J]/=E;
    }
}

void BallApproximator::calculateApproximation(const QString &resultPath, bool writeToFile)
{
    if (writeToFile)
    {
        FLIST=fopen(resultPath.toStdString().c_str(), "w");
    }

    double lenght;
    int I,J;
    double P[3], X[nonLinearParamsCount], PQ[3][maxNumberOfMeasures], timeSame[maxNumberOfMeasures];

    for (J=0; J<3; J++) DC[J]=CAM1[J]-CAM0[J];
    fisrtMeasureCount=N0; if (fisrtMeasureCount>N1) fisrtMeasureCount=N1;
    FIRST();
    if (writeToFile)
    {
        fprintf(FLIST,"\ngeometry estimatuions: t, x,y,z, error \n");
    }

    for (I=0; I<fisrtMeasureCount; I++) {
        DIST(U0[I],U1[I],P,&lenght);
        if(writeToFile){
            fprintf(FLIST,"%2d  %6.4f  %12.5e %12.5e %12.5e  %10.3e\n",I,T0[I],P[0],P[1],P[2],lenght);
        }
        //fprintf(FLIST,"%2d  %6.4f  %12.5e %12.5e %12.5e  %10.3e\n",I,T0[I],P[0],P[1],P[2],lenght);
        for (J=0; J<3; J++) PQ[J][I]=P[J];
        timeSame[I]=0.5*(T0[I]+T1[I]);
    }
    if(writeToFile){

        fprintf(FLIST,"\nlinear estimations: X[9]=(x,vx,ax/2,y,..,az/2)\n");

    }
    //fprintf(FLIST,"\nlinear estimations: X[9]=(x,vx,ax/2,y,..,az/2)\n");
    PARAB(fisrtMeasureCount,timeSame,PQ[0],P,&lenght);
    if(writeToFile){

        fprintf(FLIST,"x,vx,ax/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);

    }
    //fprintf(FLIST,"x,vx,ax/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);
    xLinearParams[0] = P[0];
    xLinearParams[1] = P[1];
    xLinearParams[2] = P[2];
    X[0]=P[0]; X[3]=P[1]; X[6]=P[2];
    PARAB(fisrtMeasureCount,timeSame,PQ[1],P,&lenght);
    if(writeToFile){

        fprintf(FLIST,"y,vy,ay/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);

    }
    //fprintf(FLIST,"y,vy,ay/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);
    yLinearParams[0] = P[0];
    yLinearParams[1] = P[1];
    yLinearParams[2] = P[2];


    X[1]=P[0]; X[4]=P[1]; X[7]=P[2];
    PARAB(fisrtMeasureCount,timeSame,PQ[2],P,&lenght);
    X[2]=P[0]; X[5]=P[1]; X[8]=P[2];
    if(writeToFile){

        fprintf(FLIST,"z,vz,az/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);

    }
    //fprintf(FLIST,"z,vz,az/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);
    zLinearParams[0] = P[0];
    zLinearParams[1] = P[1];
    zLinearParams[2] = P[2];

    if(writeToFile){

        fprintf(FLIST,"\nnonlinear estimatuions: t, x,y,z, error \n");
    }
    readData(fR, sR, fTimeR, sTimeR, fPosR, sPosR);
    //fprintf(FLIST,"\nnonlinear estimatuions: t, x,y,z, error \n");

    SOLVER(X,writeToFile);
    COVMAT(nonLinearParamsCount, 2 * N0 + 2 * N1 - nonLinearParamsCount, X,writeToFile);

    fclose(FLIST);
}



void Check(double x,double y, double z,double x1,double y1,double z1, double ax, double ay, double az, double vx, double vy, double vz,double con,
           double &t1, double &t2, double &cordx, double &cordy, double &cordz,double &cordx2,
           double &cordy2, double &cordz2) {
    double C,A,B,D;
    C = ((x1 * x) + (y1 * y) + (z1 * z)) + con;
    B = (x1 * vx) + (y1 * vy) + (z1 * vz);
    A = ((x1 * ax) ) + ((y1 * ay) ) + ((z1 * az) );
    D = sqrt((B * B) - 4 * A * C);
    if (((B * B) - (4 * A * C)) >= 0) {

        t1 = (- B + D) / (2 * A);
        t2 = (- B - D) / (2 * A);


        cordx = x + vx * (t1) + ((ax * t1 * t1));
        cordy = y + vy * (t1) + ((ay * (t1) * (t1)));
        cordz = z + vz * (t1) + ((az * (t1) * (t1)));


        cordx2 = x + vx * (t2) + ((ax * (t2) * (t2)));
        cordy2 = y + vy * (t2) + ((ay * (t2) * (t2)));
        cordz2 = z + vz * (t2) + ((az * (t2) * (t2)));}
    else {
        t1=-1;
        t2=-1;

    }

}


void calculateIntercept(double plane[4], double pos[3], double v[3], double a[3],
           double &t1, double &t2, double coord1[3], double coord2[3])
{
//    double C,A,B,D;
//    C = ((plane[0] * pos[0]) + (plane[1] * pos[1]) + (plane[2] * pos[2])) + plane[3];
//    B = (plane[0] * v[0]) + (plane[1] * v[1]) + (plane[2] * v[2]);
//    A = ((plane[0] * a[0]) ) + ((plane[1] * a[1]) ) + ((plane[2] * a[2]));
//    D = sqrt((B * B) - 4 * A * C);
//    if (((B * B) - (4 * A * C)) >= 0)
//    {

//        t1 = (- B + D) / (2 * A);
//        t2 = (- B - D) / (2 * A);


//        coord1[0] = x + vx * (t1) + ((ax * t1 * t1));
//        coord1[1] = y + vy * (t1) + ((ay * (t1) * (t1)));
//        coord1[2] = z + vz * (t1) + ((az * (t1) * (t1)));


//        coord2[0] = x + vx * (t2) + ((ax * (t2) * (t2)));
//        coord2[1] = y + vy * (t2) + ((ay * (t2) * (t2)));
//        coord2[2] = z + vz * (t2) + ((az * (t2) * (t2)));}
//    else
//    {
//        t1=-1;
//        t2=-1;

//    }

}

