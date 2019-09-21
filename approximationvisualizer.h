#ifndef APPROXIMATIONVISUALIZER_H
#define APPROXIMATIONVISUALIZER_H

#include <QObject>
#include <QGraphicsScene>
#include <calibration.h>
#include <ballapproximator.h>
#include <QImage>
#include <calibrationadjusthelper.h>
#include <QGraphicsPathItem>
#include <QPainter>
#include <qcustomplot/qcustomplot.h>
#include <qcustomplot/cxyplotter.h>
#include <QGraphicsBlurEffect>



class ApproximationVisualizer : public QObject
{
    Q_OBJECT
public:


    struct  CalibData
    {
        QDateTime dt;
        Calibration::ExteriorOr eOr;
    };

    QVector <QImage> createApproxVisualisation(BallApproximator& approx, QString info);

    void plotCamera(BallApproximator& approx, qint32 firstNumber, qint32 secondNumber, qint32 recCountF, qint32 recCountS, QDateTime dt);

    void plotCalibrate(QMap <qint32, Calibration::ExteriorOr>& map, CompareFlag compareFlag);

    void setPlots(const QVector <QCustomPlot*>& _plots)
    {
        plots = _plots;
        for (qint32 i = ObjCountGraph; i <  CalibSecondGraph + 1; ++i)
        {
            plotSyncLines.append(new QCPItemStraightLine(plots[i]));
            plotSyncLines.last()->point1->setCoords(0, 0);
            plotSyncLines.last()->point2->setCoords(0, 0);
        }
    }

    void initCamerasPlots(const QMap <qint32, QColor>& map)
    {
        colorMap = map;
        for (auto& i : colorMap.keys())
        {
            measuresMap.insert(i, QVector <Measure> ());
            calibMap.insert(i, QVector <CalibData>());
        }
    }

    void createfullPicture (bool flag) {fullPicture = flag;}

    bool isfullPicture() {return fullPicture;}

    void createshortPicture(bool flag) {shortPicture = flag;}

    bool isshortPicture() {return shortPicture;}

    qint32 getHistorySize() {return measuresMap.first().size();}

    void showMeasure(qint32 i, qint32 firstNumber, qint32 secondNumber);

    void setUpdateGraph(bool flag) {updateGraph = flag;}

    bool getUpdateGraph() {return updateGraph;}

    QImage makeFullPicture(BallApproximator& approx);

    QImage makeShortPicture(BallApproximator& approx, QString info);

    void drawTracerDebug(const QString& fVideo, const QString& sVideo);

    void clearCalibGraphs()
    {
        if (!plots.isEmpty())
        {
            plots[CalibFirstGraph]->clearGraphs();
            plots[CalibSecondGraph]->clearGraphs();
        }

    }

    static ApproximationVisualizer& instance() noexcept
    {
        static ApproximationVisualizer av;
        return av;
    }

    QMap <qint32, QVector <CalibData>> getCalibMap() {return calibMap;}

    double tBegin, tEnd, T, vBegin,
    vEnd, dxNoRot, dzNoRot, zBegin,
    xBegin, W; // tmp

    double rot[3];


signals:
    void measureCountChanged(qint32 count);
private:

    struct Measure
    {

        QVector <double> error;
        QVector <double> timeDelta;
        qint32 objCount;
        QDateTime time;
    };



    enum GraphIndex
    {
        ErrorGraph,
        TimeGraph,
        ObjCountGraph,
        CalibFirstGraph,
        CalibSecondGraph
    };

    ApproximationVisualizer(QObject* parent = nullptr);

    ~ApproximationVisualizer();


    qint32 correctCoordinates(double coord);

    void correctPixmap(QPixmap& px, double& xOld);

    void fillSkippedTime(QVector <double>& errorsFirst,  QVector <double>& timeFirstInit, QVector <double>& timeFirst);

    void clear(bool all);

    void drawLastShortPicture(bool& strike, double& vBegin, Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam,
                              Calibration::Position2D& XYShift, Calibration::Position& pShift);

    void drawInfoPlaneFullPicture(Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam,
                                  Calibration::Position2D& XYShift, Calibration::Position& pShift);

    void drawZones(Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam);



    QGraphicsScene* scenePicture = nullptr;
    QGraphicsPixmapItem* picture = nullptr;
    CXYPlotter plotter;
    bool fullPicture;
    bool shortPicture;
    bool updateGraph = true;

    QMap <qint32, QColor> colorMap;
    QVector <QCustomPlot*> plots;
    QMap <qint32, QVector <Measure>> measuresMap;
    QMap <qint32, QVector <CalibData>> calibMap;
    QVector <QCPItemStraightLine*>  plotSyncLines;
    QPixmap fullPixmap;
    QPixmap shortPixmap;

    constexpr const static qint32 ballCount = 21;
    constexpr const static double ballSize = 7.3;
    constexpr const static double coeff = 25;


    void readDrawTracerDebugData(const QString& fVideo, Calibration::ExteriorOr& EOFirstCamera,
                                 Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                                 QVector <Calibration::Position>& firstVecs, QVector <double>& firstTime, QVector<double>& videoTimes);
    void drawStrikeZoneParallelepiped(Calibration::ExteriorOr& EOFirstCamera,
                                      Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst);
};



#endif // APPROXIMATIONVISUALIZER_H
