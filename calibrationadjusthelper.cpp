#include "calibrationadjusthelper.h"

CalibrationAdjustHelper::CalibrationAdjustHelper(QObject *parent) : QObject(parent)
{

}

Mat CalibrationAdjustHelper::createCalibrateImage(const QString &refImagePath, const QString &pointListPath, Mat corImage, qint32 size)
{
    data.clear();
    Mat refImage = imread(refImagePath.toStdString(), IMREAD_GRAYSCALE);
    QFile pointListFile (pointListPath);

    if (pointListFile.open(QIODevice::ReadOnly))
    {
        QTextStream in (&pointListFile);
        QString line;
        while (in.readLineInto(&line))
        {
            auto list = line.split("\t");
            CalibratePointData d;
            d.x = list[0].toDouble();
            d.y = h - list[1].toDouble();
            d.x3d = list[2].toDouble();
            d.y3d = list[3].toDouble();
            d.z = list[4].toDouble();
            data.append(d);
        }
    }
    QVector <Mat> patterns;
    QVector <Rect> initRects;
    for (qint32 i = 0; i < data.size(); ++i)
    {
        Rect r = Rect(data[i].x - size, data[i].y - size, size * 2, size * 2);
        initRects.append(r);
        patterns.append(refImage(r));
        //imwrite(QString("/home/nvidia/recalibrate_debug/test%1.png").arg(i).toStdString(), patterns.last());
    }

    Mat corImageDraw;
    corImage.copyTo(corImageDraw);
    cvtColor(corImage, corImage, CV_BGR2GRAY);
    const qint32 searcAreaFirstCoef = 4;
    const qint32 searcAreaSecondCoef = 8;
    for (qint32 i = 0; i < patterns.size(); ++i)
    {
        rectangle(corImageDraw, initRects[i], CV_RGB(0, 255, 0));
        Rect r = Rect(data[i].x - searcAreaFirstCoef * size, data[i].y - searcAreaFirstCoef * size,
                      size * searcAreaSecondCoef, size * searcAreaSecondCoef);
        Mat img;
        corImage(r).copyTo(img);
        Mat result;
        int result_cols =  img.cols - patterns[i].cols + 1;
        int result_rows = img.rows - patterns[i].rows + 1;

        result.create( result_rows, result_cols, CV_32FC1 );

        matchTemplate( img, patterns[i], result, CV_TM_CCORR_NORMED);
        double minVal; double maxVal; Point minLoc; Point maxLoc;

        minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );

        double x = maxLoc.x + data[i].x - searcAreaFirstCoef * size;
        double y = maxLoc.y + data[i].y - searcAreaFirstCoef * size;
        rectangle( corImageDraw, Point(x, y),
                   Point( x + patterns[i].cols , y + patterns[i].rows ), CV_RGB(255, 0, 0));
        double diffX = data[i].x - size - x;
        double diffY = data[i].y - size - y;

        // qDebug() << x << data[i].x - size << y << data[i].y - size << maxLoc.x << maxLoc.y;
        data[i].x = x + size;
        data[i].y = y + size;
        putText(corImageDraw, QString("CORR%4 %1 DX: %2 DY: %3")
                .arg(maxVal)
                .arg(diffX)
                .arg(diffY)
                .arg(i).toStdString(),
                Point(x - 200, y - 10), FONT_HERSHEY_PLAIN, 1, CV_RGB(255, 255, 255));
    }
    imwrite("/home/nvidia/actual_server/server_debug/calibrate/correlate_test.png", corImageDraw);
    return corImageDraw;

}
#include <QDebug>
void CalibrationAdjustHelper::recalibrate(QVector<qint32> excludePoints, Calibration::ExteriorOr& EO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera,
                                          Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera)
{
    QVector <Calibration::Position> pos;
    QVector <Calibration::Position2D> pos2D;
    for (int i = 0; i < data.size(); ++i)
    {
        if (!excludePoints.contains(i))
        {
            Calibration::Position p;
            Calibration::Position2D pxy;
            p.X = data[i].x3d;
            p.Y = data[i].y3d;
            p.Z = data[i].z;
            pxy.X = data[i].x;
            pxy.Y = data[i].y;
            pos.append(p);
            pos2D.append(pxy);
           Calibration::Position2D d;
            Calibration::GetXYfromXYZ(EO, Camera, p, d);
          // qDebug() << "TEST" << pxy.X << d.X << pxy.Y << d.Y;
             qDebug() << pxy.X << 1216 - pxy.Y;

        }
    }
    //qDebug() << Camera.sample << Camera.line << Camera.samples << Camera.lines << Camera.focus << Camera.pixelsize << Camera.D_K1 << Camera.D_K2 << Camera.D_K3 << Camera.D_P1 << Camera.D_P2 << Camera.D_b1 << Camera.D_b2;
    testAdjustment(EO, Camera, pos, pos2D, false, newEO, newCamera);

}

void CalibrationAdjustHelper::readCurrentCalibrationParameters(qint32 cameraNumber, const QString& path, Calibration::ExteriorOr& EO,
                                                               Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera, qint32 w, qint32 h)
{
    QFile orPos (path + "/pos_orient");
    QFile otherParams (path + QString("/Cam_%1.x-cam").arg(cameraNumber));
    if (orPos.open(QIODevice::ReadOnly) && otherParams.open(QIODevice::ReadOnly))
    {
        const double pixelSize = 5.85;

        QStringList orientParams;

        QTextStream in(&orPos);
        QString line;
        while (in.readLineInto(&line))
        {
            line.replace(",", ".");
            orientParams = line.split("\t");
            if (orientParams.first().toInt() == cameraNumber)
            {
                orientParams.removeFirst();
                break;
            }
        }

        QStringList params, mainPoint;

        QXmlStreamReader xmlReader(&otherParams);
        QXmlStreamReader::TokenType type;
        while ((type = xmlReader.readNext()) != QXmlStreamReader::TokenType::EndDocument)
        {
            QXmlStreamAttributes list = xmlReader.attributes();
            if (!list.isEmpty())
            {
                if (list.hasAttribute("n"))
                {
                    if (       list.value("n") == "focus"
                               || list.value("n") == "k1"
                               || list.value("n") == "k2"
                               || list.value("n") == "k3"
                               || list.value("n") == "p1"
                               || list.value("n") == "p2"
                               || list.value("n") == "b1"
                               || list.value("n") == "b2")
                    {
                        params.append(list.value("v").toString());
                    }
                    if (list.value("n") == "principal_point")
                    {
                        for (qint32 i = 0; i < 2; ++i)
                        {
                            xmlReader.readNext();
                            mainPoint.append(xmlReader.attributes().value("v").toString());
                            xmlReader.readNext();
                        }
                    }
                }
            }
        }

        Calibration::SpacecraftPlatform::CAMERA::CameraXdirection Dir_R = Calibration::SpacecraftPlatform::CAMERA::CameraXdirection::right;
        Calibration::SpacecraftPlatform::CAMERA::CameraZdirection Dir_Z = Calibration::SpacecraftPlatform::CAMERA::CameraZdirection::frompage;
        Calibration::SpacecraftPlatform::CAMERA::CameraType CamType = Calibration::SpacecraftPlatform::CAMERA::CameraType::central;

        double pixelSizeMk = pixelSize * 1e-3;
        double mainPointX = mainPoint[0].toDouble() / pixelSizeMk  + w / 2;
        double mainPointY = mainPoint[1].toDouble() / pixelSizeMk  + h / 2;
        EO.Init(orientParams[0].toDouble(), orientParams[1].toDouble(), orientParams[2].toDouble(),
                Calibration::Deg2Radian(orientParams[3].toDouble()),
                Calibration::Deg2Radian(orientParams[4].toDouble()),
                Calibration::Deg2Radian(orientParams[5].toDouble()));
        while(params.size() < 8)
        {
            params.append("0");
        }
        Camera.Init("", params[0].toDouble(), pixelSize, mainPointX, h - mainPointY, w, h, Dir_R, Dir_Z, CamType,
                params[1].toDouble(), params[2].toDouble(), params[3].toDouble(), params[6].toDouble(), params[7].toDouble(), 0, 0,
                params[4].toDouble(), params[5].toDouble());
       // Camera.Init("", params[0].toDouble(), pixelSize, mainPointX, h - mainPointY, w, h, Dir_R, Dir_Z,CamType, Calibration::SpacecraftPlatform::CAMERA::CameraDistortionType::Classical);
    }
}

void CalibrationAdjustHelper::saveNewCalibrationParameters(Calibration::ExteriorOr &newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams &newCamera, Mat newImage, qint32 cameraId)
{

}



//newEO возвращает и принимает углы в радианах
void CalibrationAdjustHelper::testAdjustment(const Calibration::ExteriorOr& EO, const Calibration::SpacecraftPlatform::CAMERA::CameraParams& Camera,
                    QVector <Calibration::Position>& GCPs, QVector <Calibration::Position2D>& measurements,
                    bool useCameraCalibration,
                    Calibration::ExteriorOr& newEO, Calibration::SpacecraftPlatform::CAMERA::CameraParams& newCamera)
{
    //======================================================================
    Calibration::Adjustment::Data information;
    information.detectdistortion = useCameraCalibration;
    information.image_counts = 1;
    information.meas_counts = GCPs.size();
    information.images = QVector <Calibration::Adjustment::image> (information.image_counts);
    information.measurements = QVector <Calibration::Adjustment::measure> (information.meas_counts);
    //========================================================================================

    information.images[0].Camera.Init(Camera);
    information.images[0].Eo.Init(EO);
    information.images[0].ImageID = 0;
    for (int i = 0;i < information.meas_counts; i++)
    {
        information.measurements[i].camID = 0;
        information.measurements[i].GCPsID = i;
        information.measurements[i].imageID = &information.images[0];
        information.measurements[i].type = Calibration::Adjustment::pType::GCPs;
        information.measurements[i].meas = measurements[i];
        information.measurements[i].XYZg = GCPs[i];
    }
    Calibration::Adjustment ad;
    ad.Adjust(information);
    newCamera.Init(information.images[0].Camera);
    newEO.Init(information.images[0].Eo);
}

// void ExampleImprooveExtOr()
//{
//  Calibration::ExteriorOr EO = gcnew ExteriorOr();
//  Calibration::ExteriorOr newEO = gcnew ExteriorOr();
//  Calibration::SpacecraftPlatform::CAMERA::CameraParams Camera;
//  Calibration::SpacecraftPlatform::CAMERA::CameraParams newCamera;
//  //======
//  SpacecraftPlatform::CAMERA::CameraXdirection Xdir = SpacecraftPlatform::CAMERA::CameraXdirection::right;
//  SpacecraftPlatform::CAMERA::CameraZdirection Zdir = SpacecraftPlatform::CAMERA::CameraZdirection::frompage;
//  SpacecraftPlatform::CAMERA::CameraType TypeCentral = SpacecraftPlatform::CAMERA::CameraType::central;
//  SpacecraftPlatform::CAMERA::CameraDistortionType DType = SpacecraftPlatform::CAMERA::CameraDistortionType::Classical;
//  //======
//  EO.Init(-17, -5, 4.5, Deg2Radian(71.65), Deg2Radian(-67.6), Deg2Radian(-17));
//  Camera.Init("", 20.0, 5.86, 1980 / 2, 1080 / 2, 1980, 1080, Xdir, Zdir, TypeCentral, DType);
//  int counts = 9;
//  array<Position^>^ GCPs = gcnew array<Position^>(counts);
//  array<Position2D^>^ measurements = gcnew array<Position2D^>(counts);
//  int p = 0;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "A";  GCPs[p].X = 0.06;  GCPs[p].Y = 0.06;  GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Б";  GCPs[p].X = 0.41;  GCPs[p].Y = 1.12;  GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Л";  GCPs[p].X = 12.34; GCPs[p].Y = 14.11; GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "М";  GCPs[p].X = 14.11; GCPs[p].Y = 12.34; GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Н";  GCPs[p].X = 26.49; GCPs[p].Y = 1.06;  GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "П";  GCPs[p].X = 26.49; GCPs[p].Y = 26.49; GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Д";  GCPs[p].X = -1.66; GCPs[p].Y = -2.72; GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Ж";  GCPs[p].X = -2.56; GCPs[p].Y = 0.74;  GCPs[p].Z = 0.01; p++;
//  GCPs[p] = gcnew Position(); GCPs[p].name = "Г";  GCPs[p].X = 0.74;  GCPs[p].Y = -2.56; GCPs[p].Z = 0.01; p++;
//  p = 0;
//  //===== система координат связующих точек от левого верхнего угла. 1080-line - это как из системы координат PHOTOMOD перейти в СК нашей программы
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1316.86; measurements[p].Y = 1080 - 143.14; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1157.14; measurements[p].Y = 1080 - 171.43; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 334.86; measurements[p].Y = 1080 - 517.14; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 567.14; measurements[p].Y = 1080 - 528.57; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1805.14; measurements[p].Y = 1080 - 600.57; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 138.57; measurements[p].Y = 1080 - 660.57; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1775.43; measurements[p].Y = 1080 - 14.86; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1038.28; measurements[p].Y = 1080 - 37.71; p++;
//  measurements[p] = gcnew Position2D(); measurements[p].X = 1811.43; measurements[p].Y = 1080 - 129.14; p++;
//  TestAdjustment(EO, Camera, GCPs, measurements, false, newEO, newCamera);
// }
