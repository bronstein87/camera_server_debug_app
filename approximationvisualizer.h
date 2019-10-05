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
#include <limits.h>



struct PlotCameraData
{
    QVector <double> timeFirst;
    QVector <double> timeSecond;
    qint32 firstMesCount;
    qint32 secondMesCount;
};

class ApproximationVisualizer : public QObject
{
    Q_OBJECT
public:


    enum RepeatVisualizeState
    {
        ShowInit,
        ShowResultPicture,
        ShowRepeat,
        Invalid = -1
    };

    struct  CalibData
    {
        QDateTime dt;
        Calibration::ExteriorOr eOr;
    };

    static ApproximationVisualizer& instance() noexcept
    {
        static ApproximationVisualizer av;
        return av;
    }

    struct RepeatVisualizeData
    {
        bool startPointFound = false;
        bool isLastPoint = false;
        QPainterPath path;
        QPainterPath bottomPath;
        QPainterPath areaPath;
        QPainterPath linesPath;
        QRectF closeZoneEllipse;
        QRectF farZoneEllipse;
        QVector <QPolygonF> polygons;
        Calibration::Position2D p2dPrevUp;
        Calibration::Position2D p2dPrevDown;
    };



    QVector <QImage> createApproxVisualisation(BallApproximator& approx, QString info);

    void plotCamera(BallApproximator& approx, qint32 firstNumber, qint32 secondNumber, PlotCameraData& plotData, QDateTime dt);

    void plotCalibrate(QMap <qint32, Calibration::ExteriorOr>& map, CompareFlag compareFlag);

    void runRepeatStream();

    void stopRepeatStream();

    void setPlots(const QVector <QCustomPlot*>& _plots);

    void initCamerasPlots(const QMap <qint32, QColor>& map);

    void createFullPicture(bool flag) {fullPicture = flag;}

    void createShortPicture(bool flag) {shortPicture = flag;}

    bool isFullPicture() {return fullPicture;}

    bool isShortPicture() {return shortPicture;}

    qint32 getHistorySize() {return measuresMap.first().size();}

    void showMeasure(qint32 i, qint32 firstNumber, qint32 secondNumber);

    void setUpdateGraph(bool flag) {updateGraph = flag;}

    bool getUpdateGraph() {return updateGraph;}

    QImage makeFullPicture(BallApproximator& approx);

    QImage makeShortPicture(BallApproximator& approx, QString info);

    void drawTracerDebug(const QString& fVideo, const QString& sVideo);

    void clearCalibGraphs();

    double calculateBatterPositionCorr(qint32 number, cv::Mat img);

    QMap <qint32, QVector <CalibData>> getCalibMap() {return calibMap;}

    //void resetRepeatVisualise() {repeatData = RepeatVisualizeData();}

    void setCurrentRepeatState(RepeatVisualizeState state, double time = -1, cv::Mat mat = cv::Mat());

    void setLastApproximation(QSharedPointer <BallApproximator> lastApprox);

    bool appendRepeatFrame(qint32 camNum, cv::Mat frame, double time);

    void setRepeatMainInfo(qint32 camNum, double coef, double sTime);

    void setRepeatCameraNumber(qint32 number);

    void setRepeatInitTime(qint32 num, double time);

    RepeatVisualizeState getCurrentState() {return currentRepeatState;}

    bool isScenarioRun() {return repeatThreadRun;}

    void clearRepeats();

    cv::Mat drawBallTracer(RepeatVisualizeData& repeatData, QPixmap& px, Calibration::ExteriorOr& EOFirstCamera,
                           Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                           BallApproximator& approx, double videoTime, double initTime);

    double tBegin, tEnd, T, vBegin,
    vEnd, dxNoRot, dzNoRot, zBegin,
    xBegin, W, tFarZone; // tmp

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

    struct CorrTime
    {
        CorrTime() {}
        CorrTime(double _corr, double _initTime, double _startTime) :
            corr(_corr),
            initTime(_initTime),
            startTime(_startTime)
        {}
        double corr;
        double initTime;
        double startTime;
    };

    struct ScenarioData
    {
        cv::Mat initPicture;
        cv::Mat resultPicture;
        QElapsedTimer timer;
        //QLinkedList <QPair <cv::Mat, double>> repeat;
        double waitTime = -1;
        Calibration::ExteriorOr EOcamera;
        Calibration::SpacecraftPlatform::CAMERA::CameraParams camera;
        QSharedPointer <BallApproximator> approx;
        qint32 cameraNumber;
        bool corrConflictResolved = false;
        QMap <qint32, QPair <CorrTime, QLinkedList <QPair <cv::Mat, double>>>> repeats;
        QMap <qint32, cv::Rect> scaleRects;
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

    void readDrawTracerDebugData(const QString& fVideo, Calibration::ExteriorOr& EOFirstCamera,
                                 Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                                 QVector <Calibration::Position>& firstVecs, QVector <double>& firstTime, QVector<double>& videoTimes);


    void drawStrikeZoneParallelepiped(Calibration::ExteriorOr& EOFirstCamera,
                                      Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst, QGraphicsItem* parent);

    void runRepeatStreamInternal();

    void handleNewScenarioState(RepeatVisualizeState state, cv::Mat mat, double time, bool clear);



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
    QMap <qint32, cv::Rect> batterSearchZones;

    QMutex scenarioMutex;
    ScenarioData scenarioData;
    cv::VideoWriter outWindow;
    RepeatVisualizeState currentRepeatState = ShowInit;
    QAtomicInteger <qint32> repeatThreadRun = 0;

    constexpr const static qint32 ballCount = 21;
    constexpr const static double ballSize = 7.3;
    constexpr const static double coeff = 25;
    constexpr const static qint32 scenarioWidth = 1280;
    constexpr const static qint32 scenarioHeight = 720;



};



#endif // APPROXIMATIONVISUALIZER_H
