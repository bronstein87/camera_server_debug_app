#ifndef CALIBRATIONADJUSTHELPER_H
#define CALIBRATIONADJUSTHELPER_H

#include <QObject>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
//#include <opencv2/cudaimgproc.hpp>
//#include <opencv2/cudaarithm.hpp>
//#include <opencv2/cudacodec.hpp>
#include <QFile>
#include <QTextStream>
#include <calibration.h>
#include <calibrationadjusthelper.h>
#include <QXmlStreamReader>
#include <QDateTime>
#include <cameraserver.h>



struct CalibratePointData
{
    double x;
    double y;
    double x3d;
    double y3d;
    double z;
};


class CameraServer;




class CalibrationAdjustHelper : public QObject
{
    Q_OBJECT
public:

    explicit CalibrationAdjustHelper(QObject* parent = nullptr);

    cv::Mat createCalibrateImage(const QString& refImagePath, qint32 cameraNum, cv::Mat corImage, qint32 size, qint32 w, qint32 h, bool current = false);

    void recalibrate(QVector <qint32> exludePoints, Calibration::ExteriorOr& EO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& camera,
                     Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera, bool save);

    void autoCalibrate(QTcpSocket* cam, qint32 number, CameraServer* server, Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera);

    static void readCurrentCalibrationParameters(qint32 cameraNumber, const QString& path, Calibration::ExteriorOr& EO,
                                                 Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera, bool current = false,
                                                 qint32 w = 1936, qint32 h = 1216);
    void setCompareFlag(CompareFlag flag) {compareFlag = flag;}

    void setSaveAutoCalibrate(bool flag) {saveAutoCalibrate = flag;}

    void updateCalibrateInfo(qint32 cameraNum,  Calibration::ExteriorOr& eo);


private:

    void testAdjustment(const Calibration::ExteriorOr& EO, const Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera,
                        QVector <Calibration::Position>& GCPs, QVector <Calibration::Position2D>& measurements,
                        bool useCameraCalibration,
                        Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera);



    QVector <CalibratePointData> data;
    QVector <CalibratePointData> dataInit;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams newCamera;
    Calibration::ExteriorOr newEO;
    cv::Mat refImage;
    cv::Mat corImageInt;
    QVector <double> diffXVec;
    QVector <double> diffYVec;
    QVector <double> corrs;
    bool savePattens = true;
    static constexpr const qint32 w = 1936;
    static constexpr const qint32 h = 1216;
    static constexpr const qint32 whd = 1920;
    static constexpr const qint32 hhd = 1080;
    bool saveAutoCalibrate = false;
    CompareFlag compareFlag = None;

};

#endif // CALIBRATIONADJUSTHELPER_H
