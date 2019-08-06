#ifndef CALIBRATIONADJUSTHELPER_H
#define CALIBRATIONADJUSTHELPER_H

#include <QObject>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/cudaimgproc.hpp>
#include <opencv2/cudaarithm.hpp>
#include <opencv2/cudacodec.hpp>
#include <QFile>
#include <QTextStream>
#include <calibration.h>
#include <calibrationadjusthelper.h>
#include <QXmlStreamReader>


using namespace cv;
struct CalibratePointData
{
    double x;
    double y;
    double x3d;
    double y3d;
    double z;
};



class CalibrationAdjustHelper : public QObject
{
    Q_OBJECT
public:

    explicit CalibrationAdjustHelper(QObject *parent = 0);

    Mat createCalibrateImage(const QString& refImagePath, const QString& pointListPath, Mat corImage, qint32 size);

    void recalibrate(QVector <qint32> exludePoints, Calibration::ExteriorOr& EO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera,
                     Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera);

    static void readCurrentCalibrationParameters(qint32 cameraNumber, const QString& path, Calibration::ExteriorOr& EO,
                                                 Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera, qint32 w = 1936, qint32 h = 1216);

    void saveNewCalibrationParameters(Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera, Mat newImage, qint32 cameraId);


private:

    void testAdjustment(const Calibration::ExteriorOr& EO, const Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera,
                        QVector <Calibration::Position>& GCPs, QVector <Calibration::Position2D>& measurements,
                        bool useCameraCalibration,
                        Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera);

    QVector <CalibratePointData> data;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams newCamera;
    Calibration::ExteriorOr newEO;
    static constexpr const qint32 w = 1936;
    static constexpr const qint32 h = 1216;
};

#endif // CALIBRATIONADJUSTHELPER_H
