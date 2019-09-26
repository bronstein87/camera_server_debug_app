
#include "ballapproximator.h"
#include <math.h>
using namespace StrikeZone;
BallApproximator::BallApproximator()
{

}

BallApproximator::~BallApproximator()
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
    for (I=0; I<3; I++) //fprintf(FLIST,"%d %9.5f %9.5f\n",I,CAM0[I],CAM1[I]);

        if(writeToFile){
            fprintf(FLIST,"%d %9.5f %9.5f\n",I,CAM0[I],CAM1[I]);
            fprintf(FLIST,"\ndata from c0\n");
        }


    if(writeToFile){
        fprintf(FLIST,"\ndata from c0\n");
        fprintf(FLIST,"%d\n",N0);
    }

    for (I=0; I<N0; I++) fscanf(FL,"%lf %lf %lf",&U0[I][0],&U0[I][1],&U0[I][2]);
    for (I=0; I<N0; I++)
        if(writeToFile){
            fprintf(FLIST,"%2d %11.7f %11.7f %11.7f\n",I,U0[I][0],U0[I][1],U0[I][2]);
            fprintf(FLIST,"\ndata from c1\n");
            fprintf(FLIST,"%d\n",N1);
        }

    fprintf(FLIST,"\ndata from c1\n");
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

    CAM0[0] = fPos.X;
    CAM0[1] = fPos.Y;
    CAM0[2] = fPos.Z;

    CAM1[0] = sPos.X;
    CAM1[1] = sPos.Y;
    CAM1[2] = sPos.Z;

    N0 = 0;
    N1 = 0;
    qint32 i = 0;
    qint32 j = 0;
    while (i < f.size())
    {
        if (fTime[i] > 0)
        {
            ++N0;
            U0[j][0] = f[i].X;
            U0[j][1] = f[i].Y;
            U0[j][2] = f[i].Z;
            T0[j] = fTime[i];
            ++j;
        }
        ++i;
    }

    i = j = 0;
    while (i < s.size())
    {
        if (sTime[i] > 0)
        {
            ++N1;
            U1[j][0] = s[i].X;
            U1[j][1] = s[i].Y;
            U1[j][2] = s[i].Z;
            T1[j] = sTime[i];
            ++j;
        }
        ++i;
    }


    for (qint32 i = 0; i < N0; i++)
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

    for (qint32 i = 0; i < N1; i++)
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
    //    for (int i = 0; i < N0; ++i)
    //    {
    //        qDebug() << U0[i][0] << U0[i][1] << U0[i][2] << QString::number(T0[i],'g', 11);
    //    }
    //    for (int j = 0; j < N1; j++)
    //    {
    //        qDebug() << U1[j][0] << U1[j][1] << U1[j][2] << QString::number(T1[j],'g', 11);
    //    }

    TIN=T0[0];
    if (T1[0]>TIN) TIN=T1[0];
    //fprintf(FLIST,"TIN= %11.7f\n",TIN);
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


void BallApproximator::calculateIntercept(double plane[4], double pos[3], double v[3], double a[3],
double &t1, double &t2, double coord1[3], double coord2[3])
{
    double C,A,B,D;
    C = ((plane[0] * pos[0]) + (plane[1] * pos[1]) + (plane[2] * pos[2])) + plane[3];
    B = (plane[0] * v[0]) + (plane[1] * v[1]) + (plane[2] * v[2]);
    A = ((plane[0] * a[0]) ) + ((plane[1] * a[1]) ) + ((plane[2] * a[2]));
    D = sqrt((B * B) - 4 * A * C);
    if (((B * B) - (4 * A * C)) >= 0)
    {

        t1 = (- B + D) / (2 * A);
        t2 = (- B - D) / (2 * A);


        coord1[0] = pos[0] + v[0] * (t1) + ((a[0] * t1 * t1));
        coord1[1] = pos[1] + v[1] * (t1) + ((a[1] * (t1) * (t1)));
        coord1[2] = pos[2] + v[2] * (t1) + ((a[2] * (t1) * (t1)));


        coord2[0] = pos[0] + v[0] * (t2) + ((a[0] * (t2) * (t2)));
        coord2[1] = pos[1] + v[1] * (t2) + ((a[1] * (t2) * (t2)));
        coord2[2] = pos[2] + v[2] * (t2) + ((a[2] * (t2) * (t2)));
    }
    else
    {
        t1 = -1;
        t2 = -1;

    }

}

bool BallApproximator::calculatePhysicsParameters(double& tBegin, double& tEnd, double& T, double& vBegin,
                                                  double& vEnd, double& dxNoRot, double& dzNoRot, double& zBegin,
                                                  double& xBegin, double rot[3], double& W, double& tFarZone)
{
    double pos[3], v[3], a[3];
    rotateMovementParameters(pos, v , a);
    for (qint32 i = 0; i < 3; ++i)
    {
        a[i] *= 2;
    }
    double t1, t2;
    double yBegin = 16.55;
    bool strike = false;
    double widthInit = 0;
//    const double closeZoneY = 0;
//    const double farZoneY = 0.43;
//    const double width = 0.215;
//    const double minHeight = 0.473;
//    const double maxHeight= 1.045;//1.075;
    const double offset = 0.03;

    solveQuadratic(a[1] / 2, v[1], pos[1] - yBegin, t1, t2);

    xBegin = pos[0] + v[0] * t2 + a[0] / 2 * t2 * t2;
    zBegin = pos[2] + v[2] * t2 + a[2] / 2 * t2 * t2;
    double vyBegin = v[1] + a[1] * t2;
    double vxBegin = v[0] + a[0] * t2;
    double vzBegin = v[2] + a[2] * t2;
    tBegin = t2;
    vBegin = sqrt(vxBegin * vxBegin + vyBegin * vyBegin + vzBegin * vzBegin);

    //double yEnd = 0;
    solveQuadratic(a[1] / 2, v[1], pos[1] - closeZoneY, t1, t2);
    double xEnd = pos[0] + v[0] * t2 + a[0] / 2 * t2 * t2;
    double zEnd = pos[2] + v[2] * t2 + a[2] / 2 * t2 * t2;
    double vyEnd = v[1] + a[1]  * t2;
    double vxEnd = v[0] + a[0]  * t2;
    double vzEnd = v[2] + a[2]  * t2;
    tEnd = t2;
    vEnd = sqrt(vxEnd * vxEnd + vyEnd * vyEnd + vzEnd * vzEnd);

    if (xEnd <= width + widthInit + offset && xEnd >= widthInit - width - offset
            && zEnd >= minHeight - offset && zEnd <= maxHeight + offset)
    {
        strike = true;
    }

    solveQuadratic(a[1] / 2, v[1], pos[1] - farZoneY, t1, t2);
    double xFarEnd = pos[0] + v[0] * t2 + a[0] / 2 * t2 * t2;
    double zFarEnd = pos[2] + v[2] * t2 + a[2] / 2 * t2 * t2;
    tFarZone = t2;
    if (xFarEnd <= width + widthInit + offset && xFarEnd >= widthInit - width - offset
            && zFarEnd >= minHeight - offset && zFarEnd <= maxHeight + offset)
    {
        strike = true;
    }

    double vxMag = (vxBegin + vxEnd) / 2;
    double vyMag = (vyBegin + vyEnd) / 2;
    double vzMag = (vzBegin + vzEnd) / 2;
    double vMag = (vBegin + vEnd) / 2;
    T = tEnd - tBegin;

    const double g = 9.8066;
    double axMag = a[0] - a[1] * vxMag/vyMag;
    double azMag = a[2] + g - a[1] * vzMag / vyMag;
    //double deltaAx = a[0] - axMag;
    //double deltaAy = a[1] + axMag * (vxMag/vyMag) + azMag * (vzMag / vyMag);
    //double deltaAz = a[2] + g - azMag;

    double axNoRot = a[0] - axMag;
    double azNoRot = a[2] - azMag;

    double xNoRot = xBegin + vxBegin * (tEnd - tBegin) + axNoRot * pow(tEnd - tBegin, 2) / 2;
    double zNoRot = zBegin + vzBegin * (tEnd - tBegin) + azNoRot * pow(tEnd - tBegin, 2) / 2;

    dxNoRot = xEnd - xNoRot;
    dzNoRot = zEnd - zNoRot;

    double r = 7.3 / 100 / 2.0;
    double A = 3.14 * r * r;
    double rho = 1.225;
    double m = 0.143;
    double K = 1.0 / 2.0 * rho * A / m;
    //double CD = sqrt(deltaAz * deltaAz  + deltaAy * deltaAy + deltaAz * deltaAz) / (K * vyMag * vMag);
    double CL = sqrt(axMag * axMag + azMag * azMag) / (K * vMag * vMag);
    double S = -(583.0 * CL) / (2333.0 * CL - 1120.0);
    double metersToMiles =  3.6 / 1.609;
    W = (S / 8.56e-3) * vMag * metersToMiles;
    for (qint32 i = 0; i < 3; ++i)
    {
        a[i] /= 2;
    }

    double mInit[3][3] {{1, 0 , 0}, {0, 1, 0}, {0, 0, 1}};
    double mRot[3][3];
    double pNoRot[3] {xNoRot, 0, zNoRot};
    BOKZMath::rotateOZ(-45.0 * BOKZMath::degreesToRad, mInit, mRot);
    multMatrixVector(mRot, pNoRot, rot);
    return strike;
}

void BallApproximator::getPointAt(double time, double point[])
{
    point[0] = xNonLinearParams[0] + xNonLinearParams[1] * time + xNonLinearParams[2] * time * time;
    point[1] = yNonLinearParams[0] + yNonLinearParams[1] * time + yNonLinearParams[2] * time * time;
    point[2] = zNonLinearParams[0] + zNonLinearParams[1] * time + zNonLinearParams[2] * time * time;
}

void BallApproximator::solveQuadratic(double a, double b, double c, double& x1, double& x2)
{
    double discriminant, realPart, imaginaryPart;
    discriminant = b*b - 4*a*c;

    if (discriminant > 0) {
        x1 = (-b + sqrt(discriminant)) / (2*a);
        x2 = (-b - sqrt(discriminant)) / (2*a);
    }

    else if (discriminant == 0) {
        x1 = (-b + sqrt(discriminant)) / (2*a);
    }
    else {
        realPart = -b/(2*a);
        imaginaryPart =sqrt(-discriminant)/(2*a);
    }
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




void BallApproximator::calculateApproximation(const QString &resultPath, bool writeToFile)
{
    if (writeToFile)
    {
        FLIST = fopen(resultPath.toStdString().c_str(), "w");
    }

    int I,J;
    for (I=0; I<3; I++)
        if(writeToFile){
            fprintf(FLIST,"%d %9.5f %9.5f\n",I,CAM0[I],CAM1[I]);
        }
    fprintf(FLIST,"\ndata from c0\n");

    for (I=0; I<N0; I++)
        if(writeToFile){
            fprintf(FLIST,"%2d %6.4f %11.7f %11.7f %11.7f\n",I, T0[I], U0[I][0],U0[I][1],U0[I][2]);
        }

    fprintf(FLIST,"\ndata from c1\n");
    for (I=0; I<N1; I++)
        if(writeToFile){
            fprintf(FLIST,"%2d %6.4f %11.7f %11.7f %11.7f\n",I,T1[I],U1[I][0],U1[I][1],U1[I][2]);
        }

    double lenght;

    double P[3], X[nonLinearParamsCount], PQ[3][maxNumberOfMeasures];
    FIRST();
    for (J=0; J<3; J++)
    {
        DC[J]=CAM1[J]-CAM0[J];
    }
    if (writeToFile)
    {
        fprintf(FLIST,"\ngeometry estimatuions: t, x,y,z, error \n");
    }

    for (I=0; I<fisrtMeasureCount; I++) {
        DIST(U0[I],U1[I],P,&lenght);
        if(writeToFile){
            fprintf(FLIST,"%2d  %6.4f  %12.5e %12.5e %12.5e  %10.3e\n",I,T0[I],P[0],P[1],P[2],lenght);
        }
        for (J=0; J<3; J++) PQ[J][I]=P[J];
        //timeSame[I]=0.5*(T0[I]+T1[I]);
    }
    if(writeToFile){
        fprintf(FLIST,"\nlinear estimations: X[9]=(x,vx,ax/2,y,..,az/2)\n");

    }
    //fprintf(FLIST,"\nlinear estimations: X[9]=(x,vx,ax/2,y,..,az/2)\n");
    PARAB(fisrtMeasureCount, timeSame, PQ[0],P,&lenght);
    if(writeToFile){

        fprintf(FLIST,"x,vx,ax/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);

    }
    //fprintf(FLIST,"x,vx,ax/2,rms: %12.5e %12.5e %12.5e %10.3e\n",P[0],P[1],P[2],lenght);
    xLinearParams[0] = P[0];
    xLinearParams[1] = P[1];
    xLinearParams[2] = P[2];
    X[0]=P[0]; X[3]=P[1]; X[6]=P[2];
    PARAB(fisrtMeasureCount, timeSame, PQ[1],P , &lenght);
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

    SOLVER(X, writeToFile);
    COVMAT(nonLinearParamsCount, 3 * (N0 + N1) - nonLinearParamsCount, X,writeToFile);

    //NEWT(X);
    ERRORS(X);
    PROGN(X);

    fclose(FLIST);
    FLIST = nullptr;
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



void BallApproximator::PARAB(int K, double *T, double *X, double *RVW, double *E)
{
    int I;
    double S0,S1,S2,S3,S4,Y0,Y1,Y2,A,B,D,A0,A1,A2,L0,L1,L2,L3,L4,L5;
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
}



void BallApproximator::SYST(double *X, double (*B)[nonLinearParamsCountPlus])
{
    double A[3][3],P[3],S,DT,DT2;
    int I,J,K;

    for (I=0; I<N0; I++) {
        S=0; DT=T0[I]-TIN; DT2=DT*DT;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*X[J+3]+DT2*X[J+6]-CAM0[J]; S+=P[J]*P[J];
        }
        S=sqrt(S);
        for (J=0; J<3; J++) {
            P[J]/=S; B[3*I+J][9]=U0[I][J]-P[J];
        }
        for (J=0; J<3; J++) for (K=J; K<3; K++) {
            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
        }
        for (J=0; J<3; J++) for (K=0; K<3; K++) {
            B[3*I+J][K]=A[J][K]; B[3*I+J][K+3]=DT*A[J][K];
            B[3*I+J][K+6]=DT2*A[J][K];
        }
    }

    for (I=0; I<N1; I++) {
        S=0; DT=T1[I]-TIN; DT2=DT*DT;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*X[J+3]+DT2*X[J+6]-CAM1[J]; S+=P[J]*P[J];
        }
        S=sqrt(S);
        for (J=0; J<3; J++) {
            P[J]/=S; B[3*(I+N0)+J][9]=U1[I][J]-P[J];
        }
        for (J=0; J<3; J++) for (K=J; K<3; K++) {
            A[J][K]=-P[J]*P[K]/S; if (K==J) A[J][J]+=1/S; else A[K][J]=A[J][K];
        }
        for (J=0; J<3; J++) for (K=0; K<3; K++) {
            B[3*(N0+I)+J][K]=A[J][K]; B[3*(N0+I)+J][K+3]=DT*A[J][K];
            B[3*(N0+I)+J][K+6]=DT2*A[J][K];
        }
    }

}

void BallApproximator::NEWT(double *X)
{
    double POR=0,EPSI=1e-10,W,NEV,
            Q[nonLinearParamsCount],
            E[nonLinearParamsCount],
            AX[3*maxNumberOfMeasures][nonLinearParamsCountPlus];
    int NN,I,J;

    NN=0;
    do {
        NN=NN+1;
        SYST(X,AX);
        NEV=0;
        for (I=0; I<3*(N0+N1); I++) NEV+=AX[I][nonLinearParamsCount]*AX[I][nonLinearParamsCount];
        NEV=sqrt(NEV);
        MINFIT(3*(N0+N1),nonLinearParamsCount,1,1e-16,1e-30,AX,Q);
        //-----------------------------------------------------------------
        printf("singular values\n");
        fprintf(FLIST,"singular values\n");
        for (I=0; I<nonLinearParamsCount; I++) printf("%14.7e\n",Q[I]);
        for (I=0; I<nonLinearParamsCount; I++) fprintf(FLIST,"%14.7e\n",Q[I]);
        fprintf(FLIST,"\n");
        // -----------------------------------------------------------------
        for (I=0; I<nonLinearParamsCount; I++) if (Q[I]>POR) E[I]=AX[I][nonLinearParamsCount]/Q[I]; else E[I]=0;
        for (I=0; I<nonLinearParamsCount; I++) {
            Q[I]=0; for (J=0; J<3; J++) Q[I]+=AX[I][J]*E[J];
        }
        W=0;
        for (I=0; I<nonLinearParamsCount; I++) {W+=Q[I]*Q[I]; X[I]+=Q[I];}
        W=sqrt(W);
        printf("N=%2d step=%14.7e nev=%14.7e\n",NN,W,NEV);
        fprintf(FLIST,"N=%2d step=%14.7e nev=%14.7e\n",NN,W,NEV);
        for (I=0; I<nonLinearParamsCount; I++) printf("X[%d]=%17.10e;\n",I,X[I]);
        printf("\n");
        for (I=0; I<nonLinearParamsCount; I++) fprintf(FLIST,"X[%d]=%17.10e;\n",I,X[I]);
        fprintf(FLIST,"\n");
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

void BallApproximator::ERRORS(double *X)
{
    double W,Q[nonLinearParamsCount],/*E[nonLinearParamsCount],*/
            A[nonLinearParamsCount][nonLinearParamsCount],
            //B[nonLinearParamsCount][nonLinearParamsCount],
            AX[3*maxNumberOfMeasures][nonLinearParamsCountPlus];
    double P[3],S,R,DT,DT2;
    int I,J,K;

    SYST(X,AX);
    W=0;
    for (I=0; I<3*(N0+N1); I++) W+=AX[I][nonLinearParamsCount]*AX[I][nonLinearParamsCount];
    fprintf(FLIST,"FCRIT %14.7e",W);
    W=sqrt(W/(3*(N0+N1)-nonLinearParamsCount));
    printf("SIGMA(rad.) %12.5e\n",W);
    fprintf(FLIST,"  SIGMA(rad.) %14.7e\n\n",W);
    errorsFirst.clear();
    fprintf(FLIST,"residues in ort of c0: du[0..2],   distance,     estimate dr\n");
    for (I=0; I<N0; I++) {
        fprintf(FLIST,"%10.7f %10.7f %10.7f",AX[3*I][nonLinearParamsCount],AX[3*I+1][nonLinearParamsCount],AX[3*I+2][nonLinearParamsCount]);

        R=0; DT=T0[I]-TIN; DT2=DT*DT;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*X[J+3]+DT2*X[J+6]-CAM0[J]; R+=P[J]*P[J];
        }
        R=sqrt(R);
        S=0; for (J=0; J<3; J++) S+=AX[3*I+J][nonLinearParamsCount]*AX[3*I+J][nonLinearParamsCount];
        S=sqrt(S);
        fprintf(FLIST,"   %10.7f   %10.7f\n",R,R*S);
        errorsFirst.append(R * S);
    }

    fprintf(FLIST,"\n");
    errorsSecond.clear();
    fprintf(FLIST,"residues in ort of c1: du[0..2],   distance,     estimate dr\n");
    for (I=0; I<N1; I++) {
        fprintf(FLIST,"%10.7f %10.7f %10.7f",AX[3*(I+N0)][nonLinearParamsCount],AX[3*(I+N0)+1][nonLinearParamsCount],AX[3*(I+N0)+2][nonLinearParamsCount]);

        R=0; DT=T1[I]-TIN; DT2=DT*DT;
        for (J=0; J<3; J++) {
            P[J]=X[J]+DT*X[J+3]+DT2*X[J+6]-CAM1[J]; R+=P[J]*P[J];
        }
        R=sqrt(R);
        S=0; for (J=0; J<3; J++) S+=AX[3*(I+N0)+J][nonLinearParamsCount]*AX[3*(I+N0)+J][nonLinearParamsCount];
        S=sqrt(S);
        fprintf(FLIST,"   %10.7f   %10.7f\n",R,R*S);
        errorsSecond.append(R * S);
    }

    fprintf(FLIST,"\n");

    MINFIT(3*(N0+N1),nonLinearParamsCount,1,1e-16,1e-30,AX,Q);
    fprintf(FLIST,"\nsingular values and appropriate columns of V\n");
    for (I=0; I<6; I++) fprintf(FLIST,"%12.5e ",Q[I]);
    fprintf(FLIST,"\n\n");
    for (I=0; I<nonLinearParamsCount; I++) {
        for (J=0; J<6; J++) fprintf(FLIST,"%12.5e ",AX[I][J]);
        fprintf(FLIST,"\n");
    }
    fprintf(FLIST,"\n");
    for (I=6; I<nonLinearParamsCount; I++) fprintf(FLIST,"%12.5e ",Q[I]);
    fprintf(FLIST,"\n\n");
    for (I=0; I<nonLinearParamsCount; I++) {
        for (J=6; J<nonLinearParamsCount; J++) fprintf(FLIST,"%12.5e ",AX[I][J]);
        fprintf(FLIST,"\n");
    }
    fprintf(FLIST,"\n");

    for (I=0; I<nonLinearParamsCount; I++) for (J=0; J<nonLinearParamsCount; J++) {
        A[I][J]=0; for (K=0; K<nonLinearParamsCount; K++) A[I][J]+=AX[I][K]*AX[J][K]/Q[K]/Q[K];
    }
    fprintf(FLIST,"standart deviations\n");
    for (I=0; I<6; I++) fprintf(FLIST,"%12.5e ",W*sqrt(A[I][I]));
    fprintf(FLIST,"\n");
    for (I=6; I<nonLinearParamsCount; I++) fprintf(FLIST,"%12.5e ",W*sqrt(A[I][I]));
    fprintf(FLIST,"\n");
}





void BallApproximator::PROGN(double *X)
{
    int I;
    double A0,A1,A2,D,S,S1,S2;
    double E0[3],E1[3],E2[3],PZ[3];

    E0[0]=1/sqrt(2); E0[1]=E0[0]; E0[2]=0;
    E1[0]=-E0[0]; E1[1]=E0[0]; E0[2]=0;
    E2[0]=0; E2[1]=0; E2[0]=1;
    PZ[0]=0; PZ[1]=0; PZ[2]=1.5;

    fprintf(FLIST,"\nprognosis\n");
    A0=0; A1=0; A2=0;
    for (I=0; I<3; I++) {
        A0+=E0[I]*X[I+6]; A1+=E0[I]*X[I+3]; A2+=E0[I]*X[I];
    }
    fprintf(FLIST,"coefficients %14.7e %14.7e %14.7e\n",A0,A1,A2);
    fprintf(FLIST,"root ");
    D=A1*A1-4*A0*A2;
    if (fabs(A0)>1e-10) {
        S=-0.5*(A1+sqrt(D))/A0; S1=-A2/A1; S2=S1*(1-A0*S1/A1);
        fprintf(FLIST,"%14.7e %14.7e %14.7e\n",S,S1,S2);
        fprintf(FLIST,"errors %14.7e %14.7e\n",S1-S,S2-S);
    }
    else {
        S1=-A2/A1; S2=S1*(1-A0*S1/A1);
        fprintf(FLIST,"%14.7e %14.7e\n",S1,S2);
        fprintf(FLIST,"%14.7e %14.7e\n",S1-S2);
        S=S2;
    }
    S1=0; S2=0;
    for (I=0; I<3; I++) {
        D=X[I]+S*(X[I+3]+S*X[I+6]); S1+=E1[I]*D; S2+=E2[I]*D;
    }
    fprintf(FLIST,"coordinates on the target %14.7e %14.7e\n",S1,S2);

    {
        int K,NR0=16,NR=NR0*NR0,M=100; //test
        double M1,M2;

        M1=0; M2=0;
        for (I=1; I<M; I++) {
            S=0; for (K=0; K<NR; K++) {D=rand(); S+=D/32767.0-0.5;}
            S=S*sqrt(12)/NR0; //fprintf(FLIST,"%g\n",S);
            M1+=S; M2+=S*S;
        }
        M1/=M; M2=sqrt(M2/M-M1*M1);
        fprintf(FLIST,"check (0, 1)? %g %g\n",M1,M2);
    }
}




void BallApproximator::FIRST()
{
    double TB,TF,Y[maxNumberOfMeasures],P[3],E;
    int I,J;

    TB=T0[0];
    if (TB<T1[0])
    {
        TB=T1[0];
    }
    TF = T0[N0-1];

    if(TF > T1[N1-1])
    {
        TF = T1[N1-1];
    }
    E=(TF-TB)/(fisrtMeasureCount - 1);
    for (I=0; I<fisrtMeasureCount; I++)
    {
        timeSame[I]=TB+E*I;
    }

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




void BallApproximator::MINFIT(int M, int N, int P, double EPS, double TOL, double (*AB)[nonLinearParamsCountPlus], double *Q)

/*
   Computation of the matrices diag(Q), V and C such that for given real M*N
   matrix A and M*P matrix B

                   U'*A*V=diag(Q) and  U'*B=C

   with orthogonal matrices U and V. The singular values and the matrices V
   and C may be used to determine X# minimizing

                (1) normF(A*X-B) and (2) normF(X)

   with the solution

              X# = V*(Psevdo-inverse of diag(Q))*C.

   The procedure can also be used to determine the complete solution
   of underdetermined linear system, i.e. rank(A)=M<N. The array Q[0..N-1]
   represents the matrix diag(Q). A and B together are to be given as
   the first M rows of the array AB[0..M-1,0..N+P-1]. V is returned in the
   first N rows and columns of AB while C in the last P colomns of AB (if P>0).
*/
{
    int I,J,K,L,L1,NP;
    double C,F,G,H,S,X,Y,Z,E[nonLinearParamsCount];

    // Householder's reduction to bidiagonal form
    G=0.0; X=0.0; NP=N+P;
    for (I=0; I<N; I++) {
        E[I]=G; S=0.0; L=I+1;
        for (J=I; J<M; J++) S+=AB[J][I]*AB[J][I];
        if (S<TOL) G=0.0; else {
            F=AB[I][I];
            if (F<0) G=sqrt(S); else G=-sqrt(S);
            H=F*G-S; AB[I][I]=F-G;
            for (J=L; J<NP; J++) {
                S=0.0;
                for (K=I; K<M; K++) S+=AB[K][I]*AB[K][J];
                F=S/H;
                for (K=I; K<M; K++) AB[K][J]+=F*AB[K][I];
            } // end(J)
        } // end{S};
        Q[I]=G; S=0.0;
        if (I<M) for (J=L; J<N; J++) S+=AB[I][J]*AB[I][J];
        if (S<TOL) G=0.0; else {
            F=AB[I][I+1];
            if (F<0) G=sqrt(S); else G=-sqrt(S);
            H=F*G-S; AB[I][I+1]=F-G;
            for (J=L; J<N; J++) E[J]=AB[I][J]/H;
            for (J=L; J<M; J++) {
                S=0.0;
                for (K=L; K<N; K++) S+=AB[J][K]*AB[I][K];
                for (K=L; K<N; K++) AB[J][K]+=S*E[K];
            } // end{J}
        } // end{S}
        Y=fabs(Q[I])+fabs(E[I]); if (Y>X) X=Y;
    } // end{I}

    // accumulation of right-hand transformations
    for (I=N-1; I>=0; I--) {
        if (G != 0) {
            H=AB[I][I+1]*G;
            for (J=L; J<N; J++) AB[J][I]=AB[I][J]/H;
            for (J=L; J<N; J++) {
                S=0.0;
                for (K=L; K<N; K++) S+=AB[I][K]*AB[K][J];
                for (K=L; K<N; K++) AB[K][J]+=S*AB[K][I];
            } // end{J}
        } // end{G}
        for (J=L; J<N; J++) {AB[I][J]=0.0; AB[J][I]=0.0;}
        AB[I][I]=1.0; G=E[I]; L=I;
    } // end{I}
    EPS=EPS*X; //N1=N+1;
    for (I=M; I<N; I++) for (J=N; J<NP; J++) AB[I][J]=0;

    // diagonalization of the bidiogonal form
    for (K=N-1; K>=0; K--) {
        // test F splitting}
M1:   for (L=K; L>=0; L--) {
            if (fabs(E[L])<=EPS) goto M3;
            if (fabs(Q[L-1])<=EPS) goto M2;
        } // end{L}

        // cancellation of E[L] if L>1
        // cancellation
M2: C=0.0;
        S=1.0; L1=L-1;
        for (I=L; I<=K; I++) {
            F=S*E[I]; E[I]*=C;
            if (fabs(F)<=EPS) goto M3;
            G=Q[I]; H=sqrt(F*F+G*G); Q[I]=H; C=G/H; S=-F/H;
            for (J=N; J<NP; J++) {
                Y=AB[L1][J]; Z=AB[I][J];
                AB[L1][J]=Y*C+Z*S; AB[I][J]=-Y*S+Z*C;
            } // end{J}
        } //e nd{I}

        // test f convergence
M3: Z=Q[K];
        if (L==K) goto M4;

        // shift from bottom 2*2 minor
        X=Q[L]; Y=Q[K-1]; G=E[K-1]; H=E[K];
        F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2*H*Y); G=sqrt(F*F+1);
        if (F<0) F=F-G; else F=F+G;
        F=((X-Z)*(X+Z)+H*(Y/F-H))/X;

        // next QR transformation
        C=1.0; S=1.0;
        for (I=L+1; I<=K; I++) {
            G=E[I]; Y=Q[I]; H=S*G; G=C*G;
            Z=sqrt(F*F+H*H); C=F/Z; S=H/Z; E[I-1]=Z;
            F=X*C+G*S; G=-X*S+G*C; H=Y*S; Y=Y*C;
            for (J=0; J<N; J++) {
                X=AB[J][I-1]; Z=AB[J][I];
                AB[J][I-1]=X*C+Z*S; AB[J][I]=-X*S+Z*C;
            } // end{J}
            Z=sqrt(F*F+H*H); C=F/Z; S=H/Z; Q[I-1]=Z;
            F=C*G+S*Y; X=-S*G+C*Y;
            for (J=N; J<NP; J++) {
                Y=AB[I-1][J]; Z=AB[I][J];
                AB[I-1][J]=Y*C+Z*S; AB[I][J]=-Y*S+Z*C;
            } // end{J}
        } // end{I};
        E[L]=0.0; E[K]=F; Q[K]=X; goto M1;

        // convergence
M4: if (Z<0) {
            // Q[K] is made non-negative
            Q[K]=-Z;
            for (J=0; J<N; J++) AB[J][K]=-AB[J][K];
        } // end{Z}
    } //end{K}
} // end{MINFIT};


