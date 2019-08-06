#ifndef BALLAPPROXIMATOR_H
#define BALLAPPROXIMATOR_H

#include <QObject>
#include <QString>
#include <QtGlobal>
#include <calibration.h>
#include <QVector>
#include <mathfunc.h>


using namespace BOKZMath;
class BallApproximator : QObject
{
    Q_OBJECT
public:
    BallApproximator();

    void calculateApproximation(const QString& resultPath, bool writeToFile);

    void readInitData(const QString &path, double deltaT,bool writeToFile);

    void readData(QVector<Calibration::Position> f, QVector<Calibration::Position> s,
                  QVector <double> fTime, QVector <double> sTime,
                  Calibration::Position fPos, Calibration::Position sPos);


    void rotateMovementParameters(double pos[3], double v[3], double a[3]);

    void calculateIntercept(double plane[4], double pos[3], double v[3], double a[3],
               double &t1, double &t2, double coord1[3], double coord2[3]) ;

    void  getXLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = xLinearParams[i];
        }
    }



    void  getYLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = yLinearParams[i];
        }
    }


    void  getZLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = zLinearParams[i];
        }
    }

    void  getXNonLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = xNonLinearParams[i];
        }
    }


    void  getYNonLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = yNonLinearParams[i];
        }
    }
    void  getZNonLinearParameters(double res[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            res[i] = zNonLinearParams[i];
        }
    }

    double getTIN() {return TIN;}
private:

    static const int numberResults = 1;
    static const int nonLinearParamsCount = 9;
    static const int maxNumberOfMeasures = 200;
     int fisrtMeasureCount = 15;

    QVector <Calibration::Position> fR;
    QVector<Calibration::Position> sR;
    QVector <double> fTimeR;
    QVector <double> sTimeR;
    Calibration::Position fPosR;
    Calibration::Position sPosR;

    double CAM0[3], CAM1[3], DC[3];
    double T0[maxNumberOfMeasures];
    double T1[maxNumberOfMeasures];
    double U0[maxNumberOfMeasures][3];
    double U1[maxNumberOfMeasures][3];
    double timeSame[maxNumberOfMeasures];
    double xLinearParams[3];
    double yLinearParams[3];
    double zLinearParams[3];
    double xNonLinearParams[3];
    double yNonLinearParams[3];
    double zNonLinearParams[3];
    double TIN;
    int N0, N1;
    FILE *FLIST;


    void PARAB(int K, double *T, double *X, double *RVW, double *E);

    void DIST(double *U, double *V, double *POS, double *E);

    void SOLVER(double *X,bool writeToFile);

    int CHOLDET1(int N, double (*A)[nonLinearParamsCount], double *P);

    void CHOLSOL1(int N, int R, double (*A)[nonLinearParamsCount], double *P, double (*B)[numberResults], double (*X)[numberResults]);

    int CHOLINV1(int N, double (*A)[nonLinearParamsCount]);

    int CHOLINV2(int N, double *A);

    int JACOBI(int N, int EIVEC, double (*A)[nonLinearParamsCount], double (*V)[nonLinearParamsCount], double *D);

    int CHOLDET2(int N, double *A);

    void CHOLSOL2(int N, int R, double *A, double (*B)[numberResults]);

    void NORMEQ(double *X, double (*Z)[nonLinearParamsCount], double (*V)[1], double *F);

    void COVMAT(int NPAR, int NDIV, double *X,bool writeToFile);

    void FIRST();
};

#endif // BALLAPPROXIMATOR_H
