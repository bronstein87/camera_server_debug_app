#include "cameraserver.h"
using namespace cv;
CameraServer::CameraServer(QObject *parent) : QObject(parent), commandServer(new QTcpServer(this)),
    avaliableCommands {SendCameraParams, RequestCameraVideo, RequestCameraFrame, ReadyToGetStream,
                       IsServerPrepareToGetStream, GetBaseBallCoordinates,
                       GetTestDataFromClient, StartStream, StopStream, RestartCamera, AskCurrentCameraParams,
                       SendCurrentFrame, SendSaveParameters, SendRecognizeVideo,
                       EndOfMessageSign}
{
    connect(commandServer.data(), &QTcpServer::newConnection, this, &CameraServer::handleNewConnection);
    connect(commandServer.data(), &QTcpServer::acceptError, this, [this](auto errorCode)
    {
        qDebug() << errorCode;
        commandServer->errorString();
        emit this->readyMessageFromServer(commandServer->errorString());
    });
    connect(this, &CameraServer::finishedGetData, this, [this](auto socket)
    {
        cameras[socket].streamIsActive = false;
        QByteArray data;
        data.append(StopStream);
        this->fillEndOfMessage(data);
        socket->write(data);
    });

    connect(&frameEveryTimer, &QTimer::timeout, this, [this]()
    {
        if (autoCalibrate)
        {
            syncFrame(JustLast, true);
            QTimer::singleShot(500, this, [this]()
            {
                CalibrationAdjustHelper adjustHelper;
                QMutexLocker lock(&autoCalibrateMutex);
                QMap <qint32, Calibration::ExteriorOr> map;
                Calibration::ExteriorOr newEO;
                Calibration::SpacecraftPlatform::CAMERA::CameraParams newCamera;
                adjustHelper.setSaveAutoCalibrate(updateAutoCalibrate);
                adjustHelper.setCompareFlag(compareFlag);
                for (auto& i : cameras.keys())
                {
                    adjustHelper.autoCalibrate(i, cameras[i].curParams.portSendStream, this, newEO, newCamera);
                    map.insert(cameras[i].curParams.portSendStream, newEO);
                }

                lock.unlock();
                ApproximationVisualizer& av = ApproximationVisualizer::instance();
                av.plotCalibrate(map, compareFlag);
            });


        }
    });
}

void CameraServer::handleNewConnection()
{
    calibrate = false;
    QTcpSocket* clientConnection = commandServer->nextPendingConnection();
    emit connectedCameraSocket(clientConnection);
    cameras.insert(clientConnection, CameraStatus());
    QHostAddress address = clientConnection->peerAddress();
    quint16 port = clientConnection->peerPort();
    connect(clientConnection, &QAbstractSocket::disconnected,
            this, [this, clientConnection, address, port]() {
        calibrate = false;
        cameras.remove(clientConnection);

        emit  readyMessageFromServer(QString("Lost connection with: %1 : %2\n")
                                     .arg(address.toString())
                                     .arg(port), clientConnection);
        emit disconnectedCameraSocket(clientConnection);
        clientConnection->deleteLater();

    });

    connect(clientConnection, &QTcpSocket::readyRead, this, [this, clientConnection]() {handleMessageFromCamera(clientConnection);});
    qDebug() << clientConnection->peerAddress().toString().remove("::ffff:");
    emit readyMessageFromServer(QString("New connection: %1 : %2\n")
                                .arg(address.toString())
                                .arg(port), clientConnection);

}


void CameraServer::startServer()
{
    startServerInternal(commandServer.data(), 57575);
}

void CameraServer::startServerInternal(QTcpServer* server, qint32 port)
{
    if (!server->listen(QHostAddress::Any, port)) {
        emit readyMessageFromServer(QString("Unable to start the server: %1.")
                                    .arg(server->errorString()));
    }
    QString ipAddress;
    QList<QHostAddress> ipAddressesList = QNetworkInterface::allAddresses();
    // use the first non-localhost IPv4 address
    for (int i = 0; i < ipAddressesList.size(); ++i) {
        if (ipAddressesList.at(i) != QHostAddress::LocalHost &&
                ipAddressesList.at(i).toIPv4Address()) {
            ipAddress = ipAddressesList.at(i).toString();
            break;
        }
    }
    if (ipAddress.isEmpty())
        ipAddress = QHostAddress(QHostAddress::LocalHost).toString();
    emit readyMessageFromServer(QString("The server is running on\n\nIP: %1\nport: %2\n")
                                .arg(ipAddress)
                                .arg(server->serverPort()));
}


void CameraServer::syncVideo(const QString& dir, bool compress)
{
    if (!calibrate)
    {
        for (const auto& i : cameras.keys())
        {
            emit readyMessageFromServer("Время не откалибровано!", i);
        }
        return;
    }
    QMapIterator<QTcpSocket*, CameraStatus> i(cameras);
    while (i.hasNext())
    {
        i.next();
        qDebug() << "TIME " << i.value().curParams.frameRate * i.value().curParams.videoDuration;
        qint32 fc = i.value().curParams.frameRate * i.value().curParams.videoDuration;
        QtConcurrent::run(this, &CameraServer::getVideoInternal, i.key(), fc,
                          i.value().curParams.portSendStream, compress, QString("%1/video%2_%3.avi")
                          .arg(dir)
                          .arg(QDateTime::currentDateTime().toString(Qt::ISODate))
                          .arg(i.value().curParams.portSendStream));
        QByteArray data;
        data.append(RequestCameraVideo);
        fillEndOfMessage(data);
        i.key()->write(data);
    }

}

void CameraServer::checkSync()
{
    waitFor = 0;
    if (cameras.size() == 2)
    {
        QSharedPointer <QMetaObject::Connection> conn (new QMetaObject::Connection); // save this conn in class
        *conn = connect(this, &CameraServer::syncTimeGot, this, [this, conn]()
        {
            ++waitFor;
            if (waitFor == 2)
            {
                qint32 delta = abs(cameras.first().machineSyncTime.msecsTo(cameras.last().machineSyncTime));
                //                qDebug() <<  "GOT SYNC TIME" << cameras.first().syncTime << cameras.last().syncTime
                //                          << (qint64)cameras.first().syncTime - (qint64)cameras.last().syncTime << cameras.first().machineSyncTime << cameras.last().machineSyncTime
                //                          << delta;
                if (delta < machineTimeEpsilon)
                {
                    syncTime = (qint64)cameras.first().syncTime - (qint64)cameras.last().syncTime;
                    calibrate = true;
                    for (auto& i : cameras.keys())
                    {
                        emit readyMessageFromServer(QString("Время успешно синхронизировано, дельта = %1 мс")
                                                    .arg(delta), i);
                    }
                }
                else
                {
                    for (auto& i : cameras.keys())
                    {
                        emit readyMessageFromServer("Cлишком большая разница в системном времени камер. "
                                                    "Время не откалибровано.", i);
                    }

                }

                waitFor = 0;
                disconnect(*conn);
            }
        });

        QMapIterator<QTcpSocket*, CameraStatus> i(cameras);
        while (i.hasNext())
        {
            i.next();
            QtConcurrent::run(std::bind(&CameraServer::getFrameInternal, this, i.key(),
                                        i.value().curParams.portSendStream, QString("test%2_%1.bmp")
                                        .arg(i.value().curParams.portSendStream)
                                        .arg(QTime::currentTime().toString(Qt::ISODate)), true, All, false));
            QByteArray data;
            data.append(RequestCameraFrame);
            fillEndOfMessage(data);
            i.key()->write(data);
        }
    }



}

void CameraServer::syncFrame(AppendFrameRule rule, bool dontSetStream)
{
    QMapIterator<QTcpSocket*, CameraStatus> i(cameras);
    while (i.hasNext())
    {
        i.next();
        QtConcurrent::run(std::bind(&CameraServer::getFrameInternal, this, i.key(),
                                    i.value().curParams.portSendStream, QString("test%2_%1.bmp")
                                    .arg(i.value().curParams.portSendStream)
                                    .arg(QTime::currentTime().toString("hh_mm_ss")), true, rule, dontSetStream));
        QByteArray data;
        data.append(RequestCameraFrame);
        fillEndOfMessage(data);
        i.key()->write(data);
    }
}

void CameraServer::testApproximation(const QString &filePath, BallApproximator& approx, const QString& prefix)
{

    QFile file(filePath);
    Calibration::ExteriorOr EOFirstCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraFirst;

    CalibrationAdjustHelper::readCurrentCalibrationParameters(4510, "calibrate", EOFirstCamera, cameraFirst);

    QVector <QStringList>* points;
    QVector <QStringList> points1;
    QVector <QStringList> points2;
    QVector <Calibration::Position> secondVecs;
    QVector <double> secondTime;
    QVector <Calibration::Position> firstVecs;
    QVector <double> firstTime;
    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;

        while (in.readLineInto(&line))
        {
            if (line.toInt() == 3850)
            {
                points = &points2;
                continue;
            }
            else if (line.toInt() == 4510)
            {
                points = &points1;
                continue;
            }
            (*points).append(line.split(" "));
        }
    }
    else
    {
        qDebug() << file.errorString();
    }



    Calibration::RayAndPoint rp;
    Calibration::Position2D XYpix_left;

    const double divideTime = 10000000.0;
    for (qint32 i = 0; i < points1.size(); ++i)
    {
        XYpix_left.X = points1[i][0].toDouble();
        XYpix_left.Y = points1[i][1].toDouble();
        Calibration::GetRayAndPoint(EOFirstCamera, cameraFirst, XYpix_left, rp);
        qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z << QString::number(points1[i][2].toDouble()/ divideTime, 'g', 10);
        qint64 time = points1[i][2].toLongLong();

        firstVecs.append(rp.Vect);
        if (points1[i][3].toInt())
        {
            firstTime.append((double)time / divideTime);
        }
        else
        {
            firstTime.append(-(double)time / divideTime);
        }

    }

    Calibration::ExteriorOr EOSecondCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraSecond;

    CalibrationAdjustHelper::readCurrentCalibrationParameters(3850, "calibrate", EOSecondCamera, cameraSecond);

    Calibration::Position2D XYpix_right;
    for (qint32 i = 0; i < points2.size(); ++i)
    {
        XYpix_right.X = points2[i][0].toDouble();
        XYpix_right.Y = points2[i][1].toDouble();
        Calibration::GetRayAndPoint(EOSecondCamera, cameraSecond, XYpix_right, rp);
        qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z << QString::number(points2[i][2].toDouble()/ divideTime, 'g', 10);
        qint64 time = points2[i][2].toLongLong();
        secondVecs.append(rp.Vect);
        if (points2[i][3].toInt())
        {
            secondTime.append((double)time / divideTime);
        }
        else
        {
            secondTime.append(-(double)time / divideTime);
        }

        //        else
        //        {
        //            secondVecs.append(rp.Vect);
        //            secondTime.append(-1);
        //        }
    }

    //    std::reverse(firstVecs.begin(), firstVecs.end());
    //    std::reverse(secondVecs.begin(), secondVecs.end());
    //    std::reverse(firstTime.begin(), firstTime.end());
    //    std::reverse(secondTime.begin(), secondTime.end());
    bool repeat = true;
    while (repeat)
    {
        approx.readData(firstVecs, secondVecs, firstTime, secondTime, EOFirstCamera.Point, EOSecondCamera.Point);
        approx.calculateApproximation(QString("games_%1/results/result_%2.txt")
                                      .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                      .arg(QTime::currentTime().toString("hh_mm_ss")), true);
        repeat = false;
        auto v1 = approx.getFirstErrors();
        auto v2 = approx.getSecondErrors();
        if (v1.size() < 5 || v2.size() < 5)
        {
            break;
        }
        double coef = 0.05;
        if (v1.first() > coef)
        {
            for (qint32 i = 0; i < firstTime.size(); ++i)
            {
                if (firstTime[i] > 0)
                {
                    firstTime.remove(i);
                    firstVecs.remove(i);
                    repeat = true;
                    break;
                }
            }
        }
        if (v1.last() > coef)
        {
            for (qint32 i = firstTime.size() - 1; i > 0 ; --i)
            {
                if (firstTime[i] > 0)
                {
                    firstTime.remove(i);
                    firstVecs.remove(i);
                    repeat = true;
                    break;
                }
            }
        }


        if (v2.first() > coef)
        {
            for (qint32 i = 0; i < secondTime.size(); ++i)
            {
                if (secondTime[i] > 0)
                {
                    secondTime.remove(i);
                    secondVecs.remove(i);
                    repeat = true;
                    break;
                }
            }
        }
        if (v2.last() > coef)
        {
            for (qint32 i = secondTime.size() - 1; i > 0 ; --i)
            {
                if (secondTime[i] > 0)
                {
                    secondTime.remove(i);
                    secondVecs.remove(i);
                    repeat = true;
                    break;
                }
            }
        }
    }

    double pos[3], v[3], a[3];
    approx.rotateMovementParameters(pos, v ,a);
    double xp[3];
    approx.getXNonLinearParameters(xp);
    double yp[3];
    approx.getYNonLinearParameters(yp);
    double zp[3];
    approx.getZNonLinearParameters(zp);
    QVector <Calibration::Position> posXYZ;
    double tin = approx.getTIN();
    //double deltaT = secondTime[1] - secondTime[0];

    std::reverse(firstVecs.begin(), firstVecs.end());
    std::reverse(secondVecs.begin(), secondVecs.end());
    std::reverse(firstTime.begin(), firstTime.end());
    std::reverse(secondTime.begin(), secondTime.end());


    qDebug() << filePath.indexOf(QRegExp("/(?:.(?!/))+$"));
    QString time = filePath.mid(filePath.indexOf(QRegExp("/(?:.(?!/))+$")) + 1);



    QFile f3850init(QString("/home/nvidia/work_video/res/plot/3850_%1%2init").arg(time).arg(prefix));
    f3850init.open(QIODevice::WriteOnly);
    QTextStream out (&f3850init);

    for (int i = 0; i < secondTime.size(); ++i)
    {
        secondTime[i] = secondTime[i] - tin;
    }

    for (int i = 0; i < secondTime.size(); ++i)
    {
        Calibration::Position p;
        double t2 = secondTime[i] * secondTime[i];
        p.X = xp[0] + xp[1] * (secondTime[i]) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (secondTime[i]) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (secondTime[i]) + zp[2] * t2;
        posXYZ.append(p);
        Calibration::Position2D XY;
        //Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, XY);
        Calibration::GetXYfromXYZ(EOSecondCamera, cameraSecond, p, XY);
        // qDebug() << p.X << p.Y << p.Z << XY.X << XY.Y;
        //qDebug() << "Point(" << XY.X << "," << XY.Y << "),";
        out << XY.X << "\t" << XY.Y << endl;

    }
    f3850init.close();

    QFile f3850(QString("/home/nvidia/work_video/res/plot/3850_%1%2").arg(time).arg(prefix));
    f3850.open(QIODevice::WriteOnly);
    out.setDevice(&f3850);

    //    for (int i = 0; i < 20; ++i)
    //    {
    //        secondTime.prepend(secondTime.first() + deltaT);
    //    }
    //    //    tin = secondTime.first();
    //    for (int i = 0; i < 20; ++i)
    //    {
    //        secondTime.append(secondTime.last() - deltaT);
    //    }

    for (int i = 0; i < secondTime.size(); ++i)
    {
        Calibration::Position p;
        //double t2 = secondTime[i] < 0 ? -((secondTime[i]) * (secondTime[i])) : ((secondTime[i]) * (secondTime[i]));
        double t2 = secondTime[i] * secondTime[i];
        p.X = xp[0] + xp[1] * (secondTime[i]) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (secondTime[i]) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (secondTime[i]) + zp[2] * t2;
        posXYZ.append(p);
        Calibration::Position2D XY;
        //Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, XY);
        Calibration::GetXYfromXYZ(EOSecondCamera, cameraSecond, p, XY);
        // qDebug() << p.X << p.Y << p.Z << XY.X << XY.Y;
        //qDebug() << "Point(" << XY.X << "," << XY.Y << "),";
        out << XY.X << "\t" << XY.Y << endl;

    }
    f3850.close();





    QFile f4510init(QString("/home/nvidia/work_video/res/plot/4510_%1%2init").arg(time).arg(prefix));
    f4510init.open(QIODevice::WriteOnly);
    out.setDevice(&f4510init);


    for (int i = 0; i < firstTime.size(); ++i)
    {
        firstTime[i] = firstTime[i] - tin;
    }


    for (int i = 0; i < firstTime.size(); ++i)
    {
        Calibration::Position p;
        //double t2 = firstTime[i] < 0 ? -((firstTime[i]) * (firstTime[i])) : ((firstTime[i]) * (firstTime[i]));
        double t2 = firstTime[i] * firstTime[i];
        p.X = xp[0] + xp[1] * (firstTime[i]) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (firstTime[i]) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (firstTime[i]) + zp[2] * t2;
        posXYZ.append(p);
        Calibration::Position2D XY;
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, XY);
        // Calibration::GetXYfromXYZ(EOSecondCamera, cameraSecond, p, XY);
        // qDebug() << p.X << p.Y << p.Z << XY.X << XY.Y;
        //qDebug() << "Point(" << XY.X << "," << XY.Y << "),";
        out << XY.X << "\t" << XY.Y  << endl;


    }
    f4510init.close();

    QFile f4510(QString("/home/nvidia/work_video/res/plot/4510_%1%2").arg(time).arg(prefix));
    f4510.open(QIODevice::WriteOnly);
    out.setDevice(&f4510);

    //    for (int i = 0; i < 20; ++i)
    //    {
    //        firstTime.prepend(firstTime.first() + deltaT);
    //    }
    //    //    tin = firstTime.first();
    //    for (int i = 0; i < 20; ++i)
    //    {
    //        firstTime.append(firstTime.last() - deltaT);
    //    }

    for (int i = 0; i < firstTime.size(); ++i)
    {
        Calibration::Position p;
        //double t2 = firstTime[i] < 0 ? -((firstTime[i]) * (firstTime[i])) : ((firstTime[i]) * (firstTime[i]));
        double t2 = firstTime[i] * firstTime[i];
        p.X = xp[0] + xp[1] * (firstTime[i]) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (firstTime[i]) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (firstTime[i]) + zp[2] * t2;
        posXYZ.append(p);
        Calibration::Position2D XY;
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, XY);
        // Calibration::GetXYfromXYZ(EOSecondCamera, cameraSecond, p, XY);
        // qDebug() << p.X << p.Y << p.Z << XY.X << XY.Y;
        //qDebug() << "Point(" << XY.X << "," << XY.Y << "),";
        out << XY.X << "\t" << XY.Y  << endl;


    }
}

void CameraServer::enableAutoCalibrate(bool flag, qint32 syncFrameEvery, bool save, CompareFlag cFlag)
{
    autoCalibrate = flag;
    updateAutoCalibrate = save;
    compareFlag = cFlag;
    frameEveryTimer.stop();
    frameEveryTimer.setInterval(syncFrameEvery * 1e3);
    frameEveryTimer.start();
}

void CameraServer::enableRecognitionAll(bool flag)
{
    for (auto& i : cameras.keys())
    {
        CameraOptions opt;
        opt.ballRecognizeFlag = flag;
        sendParametersToCamera(opt, i);
        if (!flag)
        {
            clearRecognizeData();
        }
    }
}

bool CameraServer::handleApproximation(BallApproximator& approx, qint32 fNum, double points1[maxNumberOfMeasures][measureDim], qint32 size1, qint32 sNum, double points2[maxNumberOfMeasures][measureDim], qint32 size2)
{
    Calibration::ExteriorOr EOFirstCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraFirst;
    Calibration::ExteriorOr EOSecondCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraSecond;
    QMutexLocker lock(&autoCalibrateMutex);
    CalibrationAdjustHelper::readCurrentCalibrationParameters(fNum, "calibrate", EOFirstCamera, cameraFirst, true);
    CalibrationAdjustHelper::readCurrentCalibrationParameters(sNum, "calibrate", EOSecondCamera, cameraSecond, true);
    lock.unlock();

    QVector <Calibration::Position> secondVecs;
    QVector <double> secondTime;
    QVector <Calibration::Position> firstVecs;
    QVector <double> firstTime;

    Calibration::RayAndPoint rp;
    Calibration::Position2D XYpix_left;

    //const double divideTime = 10000000.0;
    for (qint32 i = 0; i < size1; ++i)
    {
        XYpix_left.X = points1[i][0];
        XYpix_left.Y = points1[i][1];
        Calibration::GetRayAndPoint(EOFirstCamera, cameraFirst, XYpix_left, rp);
        //qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z << QString::number(points1[i][2].toDouble()/ divideTime, 'g', 10);
        double time = points1[i][2];
        firstVecs.append(rp.Vect);
        if (qFuzzyCompare(points1[i][3], 1))
        {
            firstTime.append(time);
        }
        else
        {
            firstTime.append(-time);
        }

    }

    Calibration::Position2D XYpix_right;
    for (qint32 i = 0; i < size2; ++i)
    {
        XYpix_right.X = points2[i][0];
        XYpix_right.Y = points2[i][1];
        Calibration::GetRayAndPoint(EOSecondCamera, cameraSecond, XYpix_right, rp);
        //qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z << QString::number(points2[i][2].toDouble(), 'g', 10);
        double time = points2[i][2];
        secondVecs.append(rp.Vect);
        if (qFuzzyCompare(points2[i][3], 1))
        {
            secondTime.append(time);
        }
        else
        {
            secondTime.append(-time);
        }
    }

    approx.readData(firstVecs, secondVecs, firstTime, secondTime, EOFirstCamera.Point, EOSecondCamera.Point);
    approx.calculateApproximation(QString("games_%1/results/result%2.txt")
                                  .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                  .arg(QTime::currentTime().toString("hh_mm_ss")), true);

    auto v1 = approx.getFirstErrors();
    //qDebug() << "error1" << v1;
    double coef = 0.05;
    bool repeat = false;
    while (v1.size() > 0 && v1.first() > coef)
    {
        firstTime.removeFirst();
        firstVecs.removeFirst();
        v1.removeFirst();
        repeat = true;
    }
    while (v1.size() > 0  && v1.last() > coef)
    {
        firstTime.removeLast();
        firstVecs.removeLast();
        v1.removeLast();
        repeat = true;
    }
    auto v2 = approx.getSecondErrors();
    //qDebug() << "error2" << v2;
    while ( v2.size() > 0  &&  v2.first() > coef)
    {
        secondTime.removeFirst();
        secondVecs.removeFirst();
        v2.removeFirst();
        repeat = true;
    }
    while ( v2.size() > 0 && v2.last() > coef)
    {
        secondTime.removeLast();
        secondVecs.removeLast();
        v2.removeLast();
        repeat = true;
    }
    if  (v1.size() <= 4 || v2.size() <= 4)
    {
        for (auto& i : cameras.keys())
        {
            emit readyMessageFromServer("После отбраковки осталось слишком мало измерений. Бросок не будет обработан", i);
        }
        return false;
    }
    if (repeat)
    {
        approx.readData(firstVecs, secondVecs, firstTime, secondTime, EOFirstCamera.Point, EOSecondCamera.Point);
        approx.calculateApproximation(QString("games_%1/results/result%2.txt")
                                      .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                      .arg(QTime::currentTime().toString("hh_mm_ss")), true);

    }

    return true;
}
void CameraServer::correctSyncTime()
{
    if (cameras.first().times.size() > 0
            && cameras.last().times.size() > 0)
    {
        //qDebug() << "TIME SIZES" << cameras.first().times.size() << cameras.last().times.size();
        quint64 timeToSearch = cameras.first().times[cameras.first().times.size() / 2] - syncTime;
        auto times = cameras.last().times;
        for (qint32 i = 0; i < times.size(); ++i)
        {
            auto diff =  (qint64)timeToSearch - (qint64)times[i];
            if (abs(diff) < syncTimeEpsilon)
            {
                //qDebug() << "NEW OLD SYNC TIME " << syncTime << syncTime + diff;
                syncTime = syncTime + diff;
            }
            //qDebug() << timeToSearch << times[i] << abs(diff) << syncTimeEpsilon;
        }
        cameras.first().times.clear();
        cameras.last().times.clear();
    }
}


void CameraServer::handleMessageFromCamera(QTcpSocket* camera)
{
    QByteArray& buffer = cameras[camera].buffer;
    bool read = true;
    QByteArray endOfMessage;
    fillEndOfMessage(endOfMessage);
    while (read)
    {
        buffer.append(camera->readAll());
        qint32 pos = buffer.indexOf(endOfMessage);
        if (avaliableCommands.contains(buffer[0]) &&
                pos != -1)
        {
            qint32 packageFeature = buffer[0];
            qint32 size = pos - 1;
            auto clearData = buffer.remove(0, 1).remove(pos - 1, endOfMessage.size());

            switch (packageFeature)
            {
            case GetTestDataFromClient:
                emit readyMessageFromServer(QString("Data from: %1 : %2\n")
                                            .arg(camera->peerAddress().toString())
                                            .arg(QString(clearData.mid(0, size))), camera);
                break;
            case GetBaseBallCoordinates:
            {
                emit readyMessageFromServer(QString("Координаты мяча получены."), camera);
                RecognizeData data;
                RecognizeVideoData vData;
                quint32 size =  maxNumberOfMeasures * measureDim * sizeof(double);
                memcpy(data.data, buffer.data(), size);
                memcpy(&vData.startTime, ((char*)data.data + size  - sizeof(QTime) - sizeof(qint32)),  sizeof(QTime));
                data.startTime = vData.startTime;
                memcpy(&vData.frameCount, ((char*)data.data + size - sizeof(qint32)),  sizeof(qint32));
                memcpy(&data.size, static_cast <char*> (buffer.data()) + size, sizeof(qint32));

                QFile savef(QString("games_%2/results/points%1_%2")
                            .arg(cameras[camera].curParams.portSendStream)
                            .arg(QDate::currentDate().toString("dd_MM_yyyy")));
                savef.open(QIODevice::Append);
                QTextStream out(&savef);
                out << cameras[camera].curParams.portSendStream << "\t" << vData.startTime.toString() << endl;
                for (qint32 i = 0; i < data.size; ++i)
                {
                    quint64 timeStamp = data.data[i][2] * 1e7;
                    quint64 timeStampOld = timeStamp;
                    correctCameraTime(timeStamp, camera);
                    data.data[i][2] = (double)timeStamp / 1e7;
                    out <<  data.data[i][0] << "\t" << data.data[i][1] << "\t" << QString::number(data.data[i][2], 'g', 10)
                            << "\t" << data.data[i][3] << "\t" << timeStamp  << "\t" << timeStampOld << endl;
                }
                out << endl << endl;
                cameras[camera].recData.append(data);
                cameras[camera].recVideoData.append(vData);
                checkRecognizeResults();
                break;
            }
            case SendRecognizeVideo:
            {

               // cameras[camera].recVideoData.append(RecognizeVideoData()); // tmp
                auto& v = cameras[camera].recVideoData.first();
                //v.startTime = QTime::currentTime();
                //v.frameCount = 180; //tmp
                for (qint32 i = 0; i < v.frameCount; ++i)
                {
                    quint64 time;
                    memcpy(&time, buffer.data() + i * sizeof(quint64), sizeof(quint64));
                    v.times.append(time);
                }

               // v.frameCount = v.times.size(); // tmp
                qDebug() << QTime::currentTime() << cameras[camera].recVideoData.size() << v.frameCount;
                QtConcurrent::run(this, &CameraServer::receiveRecognizeVideo, camera);
                break;
            }
            case IsServerPrepareToGetStream:
                break;
            case AskCurrentCameraParams:
            {
                CurrentCameraParams curParams;
                memcpy(&curParams, clearData.data(), sizeof(CurrentCameraParams));
                cameras[camera].curParams = curParams;
                emit currentCameraParamsReady(camera, curParams);
                break;
            }
            case SendCurrentFrame:
            {
                quint64 time;
                memcpy(&time, clearData.data(), sizeof(quint64));
                cameras[camera].times.append(time);
                if (cameras[camera].times.size() >= correctSyncTimeEvery
                        && cameras.size() == 2 && calibrate)
                {
                    correctSyncTime();
                }
                break;
            }

            }
            buffer.remove(0, size);
        }
        else
        {
            read = false;
        }
    }

}


void CameraServer::appendFrameToBuffer(cv::Mat frame, QTcpSocket* camera, AppendFrameRule rule)
{
    QMutexLocker lock(&streamMutex);
    if (rule == All || rule == JustQueue)
    {
        cameras[camera].frameBuffer.append(frame);
    }
    if (rule != JustQueue)
    {
        cameras[camera].lastFrame = frame;
    }

    lock.unlock();
    if (rule == All || rule == JustQueue)
    {
        emit frameFromStreamReady(camera);
    }
}

void CameraServer::getFrameInternal(QTcpSocket* camera, qint32 port, const QString& path, bool sync = false, AppendFrameRule rule, bool dontSetStream)
{
    QMutexLocker lock(&streamMutex);
    if (dontSetStream || !cameras[camera].streamIsActive)
    {
        if (!dontSetStream)
        {
            cameras[camera].streamIsActive = true;
        }

        lock.unlock();
        QScopedPointer <QTcpServer> server(new QTcpServer());
        startServerInternal(server.data(), port);
        qint32 width = cameras[camera].curParams.width;
        qint32 height = cameras[camera].curParams.height;
        if (server->waitForNewConnection(10000))
        {
            qDebug() << "new client" << server->errorString();
            auto socket = server->nextPendingConnection();

            QByteArray buffer;
            QDataStream in(&buffer, QIODevice::ReadOnly);
            in.device()->reset();
            //qDebug() << "new conn" << socket;
            qint32 i = 1;
            qint32 fullSize = width * height + sizeof(QTime) + sizeof(quint64);
            while (socket->waitForReadyRead(10000))
            {
                while (socket->bytesAvailable() > 0)
                {
                    buffer.append(socket->readAll());
                    if (buffer.size() >= fullSize)
                    {
                        in.device()->reset();
                        cv::Mat matRaw (height, width, CV_8UC1);

                        quint64 intTime;
                        in.readRawData((char*)&intTime, sizeof(quint64));

                        QTime t;
                        in.readRawData((char*)&t, sizeof(QTime));
                        if (sync)
                        {
                            //qDebug() << "INT TIME" << intTime;
                            emit syncTimeGot();
                            cameras[camera].syncTime = intTime;
                            cameras[camera].machineSyncTime = t;
                        }

                        qDebug()<< "read" << in.readRawData((char*)matRaw.data, width * height) << t << intTime;
                        buffer = buffer.remove(0, fullSize);
                        Mat mat;
                        cvtColor(matRaw, mat, CV_BayerBG2BGR);
                        appendFrameToBuffer(mat, camera, rule);
                        if (!path.isEmpty())
                        {
                            cv::imwrite(path.toStdString(), mat);
                        }
                        ++i;
                    }
                }
            }
            server->close();
            if (!dontSetStream)
            {
                cameras[camera].streamIsActive = false;
            }
        }
        else
        {
            qDebug() << "ERROR" << server->errorString();
        }

    }
}

qint32 CameraServer::correctCameraTime(quint64& intTime, QTcpSocket* camera)
{
    if (cameras.firstKey() == camera && cameras.size() > 1 && calibrate)
    {
        intTime = intTime - syncTime;
    }
    qint32 exposureCenterOffset = (cameras[camera].curParams.exposure / 2) * exposureToCameraIntTime;
    intTime += exposureCenterOffset;

    return exposureCenterOffset;
}

void CameraServer::writeRawVideoStream(QTcpSocket* camera, QByteArray& buffer, cv::Mat matRaw, cv::Mat mat, qint32& i, WriteRawVideoData &d)
{
    d.in.device()->reset();
    quint64 intTime, oldIntTime;
    d.in.readRawData((char*)&intTime, sizeof(qint64));
    oldIntTime = intTime;

    qint32 exposureCenterOffset = correctCameraTime(intTime, camera);
    QTime t;
    d.in.readRawData((char*)&t, sizeof(QTime));
    qDebug() <<  camera << "read FRAME" << d.in.readRawData((char*)matRaw.data, d.width * d.height) << t << intTime << oldIntTime << syncTime;
    buffer = buffer.remove(0, d.width * d.height + sizeof(QTime) + sizeof(quint64));
    if (cameras[camera].writeVideo)
    {
        cvtColor(matRaw, mat, CV_BayerBG2BGR);
        d.video.write(mat);
    }
    if (cameras[camera].showVideo)
    {
        cvtColor(matRaw, mat, CV_BayerBG2BGR);
        appendFrameToBuffer(mat, camera);
    }
    d.ints << intTime << "\t" << t.toString("hh:mm:ss.zzz") << "\t" << exposureCenterOffset << "\t" << syncTime << endl;
    ++i;
}

void CameraServer::receiveRtspVideo(qint32 port, qint32 frameCount, const QString& path, QTcpSocket* camera, const QVector <quint64>& times, QTime throwTime)
{
    cameras[camera].streamIsActive = true;
    QString address = camera->peerAddress().toString().remove("::ffff:");
    qputenv("GST_DEBUG", "4");
    qDebug() << "try to open receive pipeline";
    QString pipeLine = QString("rtspsrc location=rtsp://%1:%2/vd latency=0 tcp-timeout=2000000 connection-speed=100000 protocols=GST_RTSP_LOWER_TRANS_TCP "
                               "! application/x-rtp,encoding-name=H265,payload=96 ! rtph265depay ! h265parse "
                               "! avdec_h265 ! videoconvert ! appsink sync = false")
            .arg(address)
            .arg(port);
    cv::VideoCapture cap(pipeLine.toStdString(), cv::CAP_GSTREAMER);
    Mat frame;
    qint32 width = cameras[camera].curParams.width;
    qint32 height = cameras[camera].curParams.height;
    qint32 frameRate = cameras[camera].curParams.frameRate;
    VideoWriter video;
    qDebug() << width << height << frameRate;
    qDebug() << "first opened";
    bool canWriteVideo = video.open(path.toStdString(), CV_FOURCC('X','V','I','D'), frameRate, Size(width, height));
    qDebug() << "all pipelines opened";
    QFile timeFile;
    QTextStream out;
    if (!times.isEmpty())
    {
        cameras[camera].writeVideo = true;
        QString timePath = path;
        timeFile.setFileName(timePath.replace(".avi", ".txt"));
        timeFile.open(QIODevice::WriteOnly);
        out.setDevice(&timeFile);
    }
    qint32 maxAttemptsCount = 2;
    qint32 attempts = 0;
    bool opened = false;
    while (!opened && attempts != maxAttemptsCount)
    {
        if (cap.isOpened())
        {
            opened = true;
            qDebug() << "stream opened";
            qint32 i = 0;
            const qint32 allowEmptyFrames = 10;
            qint32 emptyFrames = 0;
            while (true)
            {
                cap.read(frame);
                if (!frame.empty())
                {
                    emptyFrames = 0;
                    if (cameras[camera].showVideo)
                    {
                        appendFrameToBuffer(frame, camera);
                    }
                    quint64 tsMs = -1;
                    if (cameras[camera].writeVideo && canWriteVideo)
                    {
                        video.write(frame); // tmp
                        if (timeFile.isOpen())
                        {
                            tsMs = times[i];
                            correctCameraTime(tsMs, camera);
                            out << tsMs << endl;
                        }
                    }
                    ++i;
                    if ((frameCount != -1 && i >= frameCount)
                            || !cameras[camera].streamIsActive)
                    {
                        qDebug() << "FINISHED";
                        cameras[camera].streamIsActive = false;
                        break;
                    }
                }
                else
                {
                    ++emptyFrames;
                    if (allowEmptyFrames == emptyFrames)
                    {
                        cameras[camera].streamIsActive = false;
                        break;
                    }
                    qDebug() << "empty";
                }
            }
        }
        else
        {
            ++attempts;
            QThread::msleep(100);
            cap.open(pipeLine.toStdString(), cv::CAP_GSTREAMER);
        }
    }
    if (!opened)
    {
        cameras[camera].streamIsActive = false;
    }

}

void CameraServer::getVideoInternal(QTcpSocket* camera, qint32 frameCount, qint32 port, bool compress, const QString& path = QString())
{
    try
    {
        if (!cameras[camera].streamIsActive)
        {
            if (!cameras[camera].curParams.rawFrame)
            {
                receiveRtspVideo(port, frameCount, path, camera);
            }
            else
            {
                if (path.isEmpty())
                {
                    return;
                }

                QScopedPointer <QTcpServer> server(new QTcpServer());
                startServerInternal(server.data(), port);
                qint32 width = cameras[camera].curParams.width;
                qint32 height = cameras[camera].curParams.height;
                qint32 frameRate = cameras[camera].curParams.frameRate;
                VideoWriter video;
                if (compress)
                {
                    video.open(path.toStdString(), CV_FOURCC('X','V','I','D'), frameRate, Size(width, height));
                }
                else
                {
                    video.open(path.toStdString(), CV_FOURCC('D','I','B',' '), frameRate, Size(width, height));
                }
                qDebug() << "create video writer";
                cv::Mat matRaw (height, width, CV_8UC1);
                cv::Mat mat (height, width, CV_8UC3);
                if (server->waitForNewConnection(5000))
                {
                    auto socket = server->nextPendingConnection();
                    // qDebug() << "new conn" << socket;
                    QByteArray buffer;
                    QDataStream in(&buffer, QIODevice::ReadOnly);
                    QString reportPath = path;
                    qint32 pos = reportPath.indexOf(".");
                    reportPath.remove(pos, reportPath.size() - pos);
                    QFile f (reportPath.append(".txt"));
                    qDebug() << reportPath << f.open(QIODevice::WriteOnly)  << f.errorString();
                    QTextStream outTime (&f);
                    in.device()->reset();

                    qint32 i = 1;
                    WriteRawVideoData d(in, outTime, video, width, height);

                    while (socket->waitForReadyRead(5000))
                    {
                        while (socket->bytesAvailable() > 0)
                        {
                            buffer.append(socket->readAll());
                            if (buffer.size() >= width * height + sizeof(QTime) + sizeof(quint64))
                            {
                                writeRawVideoStream(camera, buffer, matRaw, mat, i, d);
                            }
                        }
                    }
                    while (buffer.size() > 0) // write everything that left
                    {
                        writeRawVideoStream(camera, buffer, matRaw, mat, i, d);
                    }
                    qDebug() << "finished" << i << buffer.size();
                    server->close();
                }
            }
        }
    }
    catch(std::exception& e)
    {
        qDebug() << e.what();
    }

}


void CameraServer::sendParametersToCamera(const CameraOptions& parameters, QTcpSocket* socket)
{

    if (cameras.contains(socket))
    {
        QByteArray data;
        data.append(SendCameraParams);
        data.append((const char*)&parameters, sizeof(CameraOptions));
        fillEndOfMessage(data);
        socket->write(data);
    }
}

void CameraServer::receiveRecognizeVideo(QTcpSocket* camera)
{
    //qDebug() << "prepare" << cameras[camera].streamIsActive;
    if (!cameras[camera].streamIsActive)
    {
        cameras[camera].streamIsActive = true;
        auto data = cameras[camera].recVideoData.first();
        cameras[camera].recVideoData.removeFirst();
        QString path = QString("games_%1/%2/video_%3.avi")
                .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                .arg(cameras[camera].curParams.portSendStream)
                .arg(data.startTime.toString("hh_mm_ss"));
        //qDebug() << "start";
        receiveRtspVideo(cameras[camera].curParams.portSendStream + 1,
                         data.frameCount, path,
                         camera, data.times, data.startTime);
    }
}

void CameraServer::fillEndOfMessage(QByteArray& array)
{
    for (qint32 i = 0; i < 10; ++i)
    {
        array.append(EndOfMessageSign);
    }

}

void CameraServer::requestCurrentCameraParameters(QTcpSocket* socket)
{
    QByteArray array;
    array.append(AskCurrentCameraParams);
    fillEndOfMessage(array);
    socket->write(array);
}

void CameraServer::updateParamsForAllCameras()
{
    for (auto& i : cameras.keys())
    {
        requestCurrentCameraParameters(i);
    }
}


void CameraServer::requestFrameFromCamera(QTcpSocket* socket, qint32 port, const QString& path)
{
    if (cameras.contains(socket))
    {
        QtConcurrent::run(std::bind(&CameraServer::getFrameInternal, this, socket, port, path, false, All, false));
        QByteArray data;
        data.append(RequestCameraFrame);
        fillEndOfMessage(data);
        socket->write(data);
    }
}


void CameraServer::requestVideoFromCamera(qint32 port, qint32 duration, bool compress, QTcpSocket* socket, const QString& path)
{
    if (cameras.contains(socket))
    {
        QtConcurrent::run(this, &CameraServer::getVideoInternal, socket, duration * cameras[socket].curParams.frameRate, port, compress, path);
        QByteArray data;
        data.append(RequestCameraVideo);
        fillEndOfMessage(data);
        socket->write(data);
    }
}


void CameraServer::startStream(qint32 port, QTcpSocket* socket, bool start)
{
    if (cameras.contains(socket))
    {
        QByteArray data;

        if (start && !cameras[socket].streamIsActive)
        {
            data.append(StartStream);
            fillEndOfMessage(data);
            socket->write(data);
            emit streamFromIsComing(socket);
            QtConcurrent::run(this, &CameraServer::getVideoInternal, socket, -1, port, false, QString());
        }
        else
        {
            data.append(StopStream);
            fillEndOfMessage(data);
            socket->write(data);
            cameras[socket].streamIsActive = false;

        }
    }
}

void CameraServer::restartCamera(QTcpSocket* socket)
{
    QByteArray data;
    data.append(RestartCamera);
    fillEndOfMessage(data);
    socket->write(data);
}

void CameraServer::saveCameraSettings(QTcpSocket* camera)
{
    QByteArray data;
    data.append(SendSaveParameters);
    fillEndOfMessage(data);
    camera->write(data);
}

void CameraServer::checkRecognizeResults()
{
    if (cameras.size() > 1)
    {
        if (calibrate)
        {
            QVector <RecognizeData>& resFirst = cameras.first().recData;
            QVector <RecognizeData>& resSecond = cameras.last().recData;
            double min = 100;
            qint32 indexf = -1;
            qint32 indexs = -1;
            for (qint32 i = 0; i < resFirst.size(); ++i)
            {
                for (qint32 j = 0; j < resSecond.size(); ++j)
                {
                    qDebug() << i << j << QString::number(resFirst[i].data[0][2], 'g', 10)
                            << QString::number(resSecond[j].data[0][2], 'g', 10)
                            << std::abs(resFirst[i].data[0][2] - resSecond[j].data[0][2]) << "REC TIME EVALUATE";
                    double delta = std::abs(resFirst[i].data[0][2] - resSecond[j].data[0][2]);
                    if (min > delta)
                    {
                        min = delta;
                        indexf = i;
                        indexs = j;
                    }
                }
            }
            if (min < 0.5)
            {
                for (auto& cam : cameras.keys())
                {
                    emit readyMessageFromServer(QString("Измерения синхронизированы, дельта = %1 сек.").arg(min), cam);
                }
                BallApproximator approx;
                if (handleApproximation(approx, cameras.first().curParams.portSendStream, resFirst[indexf].data, resFirst[indexf].size,
                                        cameras.last().curParams.portSendStream, resSecond[indexs].data, resSecond[indexs].size))
                {
                    ApproximationVisualizer& visualizer = ApproximationVisualizer::instance();

                    qint32 recCountFirst = 0;;
                    for (qint32 i = 0; i < resFirst[indexf].size; ++i)
                    {
                        if (resFirst[indexf].data[3])
                        {
                            ++recCountFirst;
                        }
                    }
                    qint32 recCountSecond  = 0;
                    for (qint32 i = 0; i < resSecond[indexs].size; ++i)
                    {
                        if (resSecond[indexs].data[3])
                        {
                            ++recCountSecond;
                        }
                    }
//                    QString info = QString ("%1 %2 %3 %4")
//                            .arg(cameras.first().curParams.portSendStream).arg(recCountFirst)
//                            .arg(cameras.last().curParams.portSendStream).arg(recCountSecond);
                    QVector <QImage> result = visualizer.createApproxVisualisation(approx, QString());
                    double vBeginMiles = visualizer.vBegin * metersToMiles;
                    double vEndMiles = visualizer.vEnd * metersToMiles;
                    if (vBeginMiles < vEndMiles || vBeginMiles > 100 || vEndMiles > 100)
                    {
                        for (auto& cam : cameras.keys())
                        {
                            emit readyMessageFromServer(QString("Измерения отброшены по скорости. (%1)").arg(vBeginMiles), cam);
                        }
                        resFirst.remove(0, indexf + 1);
                        resSecond.remove(0, indexs + 1);
                        return;
                    }
                    emit resultPictureReady(result.first());
                    if (recCountFirst < 16)
                    {

                        QFile savef(QString("games_%2/results/points%1_%2")
                                    .arg(cameras.first().curParams.portSendStream)
                                    .arg(QDate::currentDate().toString("dd_MM_yyyy")));
                        savef.open(QIODevice::Append);
                        QTextStream out(&savef);
                        out << "ПОСМОТРЕТЬ" << endl;
                    }
                    if (recCountSecond < 16)
                    {

                        QFile savef(QString("games_%2/results/points%1_%2")
                                    .arg(cameras.last().curParams.portSendStream)
                                    .arg(QDate::currentDate().toString("dd_MM_yyyy")));
                        savef.open(QIODevice::Append);
                        QTextStream out(&savef);
                        out << "ПОСМОТРЕТЬ" << endl;
                    }
                    if (result.size() == 2)
                    {
                        result.last().save(QString("games_%1/pictures/small.jpg").arg(QDate::currentDate().toString("dd_MM_yyyy")));
                        result.first().save(QString("games_%1/pictures/big%2.jpg")
                                            .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                            .arg(QTime::currentTime().toString("hh_mm_ss")));
                        result.last().save(QString("games_%1/pictures/small%2.jpg")
                                           .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                           .arg(QTime::currentTime().toString("hh_mm_ss")));
                        result.first().save(QString("games_%1/pictures/big.jpg").arg(QDate::currentDate().toString("dd_MM_yyyy")));

                    }
                    else
                    {
                        if (visualizer.isfullPicture())
                        {
                            result.first().save(QString("games_%1/pictures/big.jpg").arg(QDate::currentDate().toString("dd_MM_yyyy")));
                            result.first().save(QString("games_%1/pictures/big%2.jpg")
                                                .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                                .arg(QTime::currentTime().toString("hh_mm_ss")));
                        }
                        else
                        {
                            result.last().save(QString("games_%1/pictures/small.jpg").arg(QDate::currentDate().toString("dd_MM_yyyy")));
                            result.first().save(QString("games_%1/pictures/small%2.jpg")
                                                .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                                .arg(QTime::currentTime().toString("hh_mm_ss")));
                        }
                    }




                    QDateTime dt = QDateTime::currentDateTime();
                    dt.setTime(resFirst[indexf].startTime);
                    visualizer.plotCamera(approx, cameras.first().curParams.portSendStream, cameras.last().curParams.portSendStream,
                                          recCountFirst, recCountSecond, dt);
                }

                resFirst.remove(0, indexf + 1);
                resSecond.remove(0, indexs + 1);
            }
            else
            {
                for (auto& cam : cameras.keys())
                {
                    emit readyMessageFromServer(QString("Не обнаружено синхронного измерения, min = %1").arg(min), cam);
                }
            }
        }

        else
        {
            emit readyMessageFromServer("Время не синхронизировано. Измерения не будут обработаны");
        }
    }
}


cv::Mat CameraServer::getFrameFromStream(QTcpSocket* socket)
{
    QMutexLocker lock(&streamMutex);
    QLinkedList <cv::Mat>& list = cameras[socket].frameBuffer;
    if (!list.isEmpty())
    {
        auto frame = list.first();
        list.pop_front();
        return frame;
    }

    return cv::Mat();
}

Mat CameraServer::getLastCameraFrame(QTcpSocket* camera)
{
    QMutexLocker lock(&streamMutex);
    return cameras[camera].lastFrame;
}

Mat CameraServer::getLastCameraFrame(qint32 num)
{
    QMutexLocker lock(&streamMutex);
    for (auto& i : cameras)
    {
        if (i.curParams.portSendStream == num)
        {
            return i.lastFrame;
        }
    }
    return Mat();
}


void CameraServer::createVideoTimer(qint32 interval, qint32 duration, QTime from, QTime to, qint32 port, QTcpSocket* socket, bool compress, const QString& path)
{
    QTimer*& timer = cameras[socket].timer;
    if (timer != nullptr)
    {
        emit readyMessageFromServer("Съемка видео по таймеру уже была активирована!");
    }
    timer = new QTimer();
    connect(timer, &QTimer::timeout, [this, port, duration, path, from, to, socket, compress]()
    {
        QTime currentTime = QTime::currentTime();
        qDebug() << from << to << currentTime << (currentTime > from && currentTime < to);
        if (currentTime > from && currentTime < to)
        {
            QString currentVideoPath = path + QString("/video_%1.avi").arg(QDateTime::currentDateTime().
                                                                           toString(Qt::ISODate));
            qDebug() << currentVideoPath;
            requestVideoFromCamera(port, duration, compress, socket, currentVideoPath);
        }

    });
    timer->start(interval * 60 * 1000);
}

void CameraServer::stopVideoTimer(QTcpSocket* socket)
{
    delete cameras[socket].timer;
}


CameraServer::~CameraServer()
{
    for (auto& i : cameras.keys())
    {
        saveCameraSettings(i);
    }
    for (auto& i : cameras)
    {
        delete i.timer;
    }
    cameras.clear();
}
