#include "cameraserver.h"

CameraServer::CameraServer(QObject *parent) : QObject(parent), commandServer(new QTcpServer(this)),
    avaliableCommands {SendCameraParams, RequestCameraVideo, RequestCameraFrame, ReadyToGetStream,
                       IsServerPrepareToGetStream, GetBaseBallCoordinates,
                       GetTestDataFromClient, StartStream, StopStream, RestartCamera, AskCurrentCameraParams, SendCurrentFrame, EndOfMessageSign}
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
        data.append(EndOfMessageSign);
        socket->write(data);
    });
}

void CameraServer::handleNewConnection()
{
    calibrate = false;
    QTcpSocket* clientConnection = commandServer->nextPendingConnection();
    emit connectedCameraSocket(clientConnection);
    qDebug() << clientConnection;
    cameras.insert(clientConnection, CameraStatus());
    QHostAddress address = clientConnection->peerAddress();
    qint16 port = clientConnection->peerPort();
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
        data.append(EndOfMessageSign);
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
                qDebug() <<  "GOT SYNC TIME" << cameras.first().syncTime << cameras.last().syncTime
                          << (qint64)cameras.first().syncTime - (qint64)cameras.last().syncTime << cameras.first().machineSyncTime << cameras.last().machineSyncTime
                          << abs(cameras.first().machineSyncTime.msecsTo(cameras.last().machineSyncTime));
                if (abs(cameras.first().machineSyncTime.msecsTo(cameras.last().machineSyncTime)) < machineTimeEpsilon)
                {
                    syncTime = (qint64)cameras.first().syncTime - (qint64)cameras.last().syncTime;
                    calibrate = true;
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
            QtConcurrent::run(this, &CameraServer::getFrameInternal, i.key(),
                              i.value().curParams.portSendStream, QString("test%2_%1.bmp")
                              .arg(i.value().curParams.portSendStream)
                              .arg(QTime::currentTime().toString(Qt::ISODate)), true);
            QByteArray data;
            data.append(RequestCameraFrame);
            data.append(EndOfMessageSign);
            i.key()->write(data);
        }
    }



}

void CameraServer::syncFrame()
{
    QMapIterator<QTcpSocket*, CameraStatus> i(cameras);
    while (i.hasNext())
    {

        i.next();

        QtConcurrent::run(this, &CameraServer::getFrameInternal, i.key(),
                          i.value().curParams.portSendStream, QString("test%2_%1.bmp")
                          .arg(i.value().curParams.portSendStream)
                          .arg(QTime::currentTime().toString(Qt::ISODate)), true);
        QByteArray data;
        data.append(RequestCameraFrame);
        data.append(EndOfMessageSign);
        i.key()->write(data);
    }
}

void CameraServer::testApproximation(const QString &fCameraPath, const QString &sCameraPath, QGraphicsPixmapItem* item)
{

    QFile fCameraPoints(fCameraPath);
    Calibration::ExteriorOr EO_4510_left;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams Camera_4510_left;

    CalibrationAdjustHelper::readCurrentCalibrationParameters(4510, "/home/nvidia/actual_server/server_debug/calibrate", EO_4510_left, Camera_4510_left);

    QStringList orientParams;
    QVector <QStringList> points;
    if (fCameraPoints.open(QIODevice::ReadOnly))
    {
        QTextStream in(&fCameraPoints);
        QString line;
        while (in.readLineInto(&line))
        {
            points.append(line.split(" "));
        }
    }


    QVector <Calibration::Position> firstVecs;
    QVector <double> firstTime;
    Calibration::RayAndPoint rp;
    Calibration::Position2D XYpix_left;

    const double divideTime = 10000000.0;
    for (qint32 i = 0; i < points.size(); ++i)
    {
        XYpix_left.X = points[i][0].toDouble();
        XYpix_left.Y = points[i][1].toDouble();
        Calibration::GetRayAndPoint(EO_4510_left, Camera_4510_left, XYpix_left, rp);
        qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z ;
        firstVecs.append(rp.Vect);
        firstTime.append(points[i][2].toDouble() / divideTime);

    }


    QFile sCameraPoints(sCameraPath);

    orientParams.clear();
    QVector <QStringList> points2;
    if (sCameraPoints.open(QIODevice::ReadOnly))
    {
        QTextStream in(&sCameraPoints);
        QString line;
        while (in.readLineInto(&line))
        {
            points2.append(line.split(" "));
        }
    }


    Calibration::ExteriorOr EO_3850_right;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams Camera_3850_right;

    CalibrationAdjustHelper::readCurrentCalibrationParameters(3850, "/home/nvidia/actual_server/server_debug/calibrate", EO_3850_right, Camera_3850_right);


    QVector <Calibration::Position> secondVecs;
    QVector <double> secondTime;
    Calibration::Position2D XYpix_right;
    for (qint32 i = 0; i < points2.size(); ++i)
    {
        XYpix_right.X = points2[i][0].toDouble();
        XYpix_right.Y = points2[i][1].toDouble();
        Calibration::GetRayAndPoint(EO_3850_right, Camera_3850_right, XYpix_right, rp);
        qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z ;
        secondVecs.append(rp.Vect);
        secondTime.append((double)points2[i][2].toDouble() / divideTime);
    }

    BallApproximator approx;
    std::reverse(firstVecs.begin(), firstVecs.end());
    std::reverse(secondVecs.begin(), secondVecs.end());
    std::reverse(firstTime.begin(), firstTime.end());
    std::reverse(secondTime.begin(), secondTime.end());
    approx.readData(firstVecs, secondVecs, firstTime, secondTime, EO_4510_left.Point, EO_3850_right.Point);
    approx.calculateApproximation("/home/nvidia/result.txt", true);
    double pos[3], v[3], a[3];
    approx.rotateMovementParameters(pos, v ,a);
    qDebug() << "finished";
    double xp[3];
    approx.getXNonLinearParameters(xp);
    double yp[3];
    approx.getYNonLinearParameters(yp);
    double zp[3];
    approx.getZNonLinearParameters(zp);
    QVector <Calibration::Position> posXYZ;
    double tin = approx.getTIN();
    qDebug() << secondTime.size();
    double deltaT = secondTime[1] - secondTime[0];

    std::reverse(firstVecs.begin(), firstVecs.end());
    std::reverse(secondVecs.begin(), secondVecs.end());
    std::reverse(firstTime.begin(), firstTime.end());
    std::reverse(secondTime.begin(), secondTime.end());

    for (int i = 0; i < secondTime.size(); ++i)
    {
        secondTime[i] = secondTime[i] - tin;
    }
    //    for (int i = 0; i < 20; ++i)
    //    {
    //        secondTime.prepend(secondTime.first() + deltaT);
    //    }
    //    //    tin = secondTime.first();
    //    for (int i = 0; i < 20; ++i)
    //    {
    //        secondTime.append(secondTime.last() - deltaT);
    //    }

    //    for (int i = 0; i < secondTime.size(); ++i)
    //    {
    //        Calibration::Position p;
    //        //double t2 = secondTime[i] < 0 ? -((secondTime[i]) * (secondTime[i])) : ((secondTime[i]) * (secondTime[i]));
    //        double t2 = secondTime[i] * secondTime[i];
    //        p.X = xp[0] + xp[1] * (secondTime[i]) + xp[2] * t2;
    //        p.Y = yp[0] + yp[1] * (secondTime[i]) + yp[2] * t2;
    //        p.Z = zp[0] + zp[1] * (secondTime[i]) + zp[2] * t2;
    //        posXYZ.append(p);
    //        Calibration::Position2D XY;
    //        //Calibration::GetXYfromXYZ(EO_4510_left, Camera_4510_left, p, XY);
    //        Calibration::GetXYfromXYZ(EO_3850_right, Camera_3850_right, p, XY);
    //        // qDebug() << p.X << p.Y << p.Z << XY.X << XY.Y;
    //        qDebug() << "Point(" << XY.X << "," << XY.Y << "),";


    //    }

    Calibration::ExteriorOr eOr;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cam;
    CalibrationAdjustHelper::readCurrentCalibrationParameters(9999, "/home/nvidia/actual_server/server_debug/calibrate", eOr, cam, 1920, 1080);
    Calibration::Position2D XY;




    QPixmap ball("/home/nvidia/actual_server/server_debug/calibrate/baseball_PNG18979.png");

    for (int i = 0; i < secondTime.size(); ++i)
    {
        Calibration::Position p;
        //double t2 = secondTime[i] < 0 ? -((secondTime[i]) * (secondTime[i])) : ((secondTime[i]) * (secondTime[i]));
        double t2 = secondTime[i] * secondTime[i];
        p.X = xp[0] + xp[1] * (secondTime[i]) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (secondTime[i]) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (secondTime[i]) + zp[2] * t2;

        Calibration::GetXYfromXYZ(eOr, cam, p, XY);
        qDebug() << XY.X << XY.Y;
        Calibration::Position2D dXY;
        dXY = XY;
        p.X += 0.073;
        Calibration::GetXYfromXYZ(eOr, cam, p, XY);
        qDebug() << XY.X << XY.Y << abs(XY.X - dXY.X);
        qint32 bW = abs(XY.X - dXY.X);
        p.X -= 0.073;
        p.Z += 0.073;
        Calibration::GetXYfromXYZ(eOr, cam, p, XY);
        qDebug() << XY.X << XY.Y << abs(XY.Y - dXY.Y);
        qint32 bH = abs(XY.Y - dXY.Y);

        QGraphicsPixmapItem* cItem = new QGraphicsPixmapItem(item);

        QPixmap tmp = ball.copy();

        cItem->setPixmap(tmp.scaled(bW * 2.5, bW * 2.5));
        cItem->setPos(dXY.X , dXY.Y );
        qDebug() << "////////////////////" << dXY.X + bW * 2 << dXY.Y + bH * 2 <<ball.isNull();
    }
    item->pixmap().save("/home/nvidia/actual_server/server_debug/calibrate/pic_tst.jpg", "JPG");




    //    Calibration::Position2D XYpix_left;
    //    Calibration::Position2D XYpix_right;
    //    Calibration::Position XYZ;
    //    //===================================================================  Set EO (xyz_opk in radians)  ===================================================================
    //    EO_left.Init(-17.991471, -5.003513, 4.498402,Calibration::Deg2Radian(71.565202039601), Calibration::Deg2Radian(-67.699490387514), Calibration::Deg2Radian(-17.151111672969));
    //    EO_right.Init(-5.011129, -17.996709, 4.503400,Calibration::Deg2Radian(82.603784380449), Calibration::Deg2Radian(-21.118703676195), Calibration::Deg2Radian(-2.667875084203));
    //    //===================================================================  Set camera   ===================================================================
    //    Camera_4510.Init("Камера 1", 2.0759815731973223e+01, 5.86, 1980 / 2, 1080 / 2, 1980, 1080, Dir_R, Dir_Z, CamType, DistT);
    //    Camera_3850.Init("Камера 2", 2.0688352029659185e+01, 5.86, 1980 / 2, 1080 / 2, 1980, 1080, Dir_R, Dir_Z, CamType, DistT);

    //    QVector<Calibration::Position2D> pos_L;
    //    QVector<baseball::Position2D> pos_R;
    //    QFile f ("C:/work/ball throw model/coords.txt");

    //        QTextStream in(&f);
    //        in.readLine();
    //        for (int i = 0; i < 16; ++i)
    //        {
    //            QString str = in.readLine();
    //            auto list = str.split("\t", QString::SkipEmptyParts);

    //            pos_L[i].X = list[2].toDouble(); pos_L[i].Y = /*1080 - */list[3].toDouble(); pos_R[i].X = list[0].toDouble(); pos_R[i].Y = /*1080 -*/ list[1].toDouble();
    //            XYpix_left = pos_L[i]; XYpix_right = pos_R[i];
    //            Calibration::RayAndPoint rp;
    //            Calibration::GetRayAndPoint(EO_left, Camera_left, XYpix_left, rp);
    //Calibration::GetRayAndPoint(EO_right, Camera_right, XYpix_right, rp);
    //           qDebug() << rp.Pos.X  << rp.Pos.Y << rp.Pos.Z << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z ;
    //            //======================================================================================================================================
    //            if (baseball::GetXYZfromXYXY(EO_left, EO_right, Camera_left, Camera_right, XYpix_left, XYpix_right, XYZ))
    //            {
    //                //baseball::RayAndPoint rp;
    //                //qDebug() << "Res: X: " + QString::number(XYZ.X) + " Y: " + QString::number(XYZ.Y) + " Z: " + QString::number(XYZ.Z);
    //                baseball::Position2D XYpix_check;
    //                if (baseball::GetXYfromXYZ(EO_left, Camera_left, XYZ, XYpix_check))
    //                {
    //                    // qDebug() <<  "Res: delta Xpix: " + QString::number(XYpix_check.X - XYpix_left.X) + " delta Ypix: " + QString::number(XYpix_check.Y - XYpix_left.Y);
    //                }
    //                if (baseball::GetXYfromXYZ(EO_right, Camera_right, XYZ, XYpix_check))
    //                {
    //                    // qDebug() << " Res: delta Xpix: " + QString::number(XYpix_check.X - XYpix_right.X) + " delta Ypix: " + QString::number(XYpix_check.Y - XYpix_right.Y);
    //                }
    //            }
    //        }
}


void CameraServer::correctSyncTime()
{
    if (cameras.size() == 2)
    {
        //qDebug() << "TIME SIZES" << cameras.first().times.size() << cameras.last().times.size();
        quint64 timeToSearch = cameras.first().times[cameras.first().times.size() / 2] - syncTime;
        auto times = cameras.last().times;
        for (qint32 i = 0; i < times.size(); ++i)
        {
            auto diff =  (qint64)timeToSearch - (qint64)times[i];
            if (abs(diff) < syncTimeEpsilon)
            {
                qDebug() << "NEW OLD SYNC TIME " << syncTime << syncTime + diff;
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
    auto buffer = cameras[camera].buffer;
    buffer.append(camera->readAll());
    //auto hexBuffer = buffer.toHex();

    if (avaliableCommands.contains(buffer[0]) &&
            buffer[buffer.size() - 1] == EndOfMessageSign)
    {
        qint32 packageFeature = buffer[0];
        auto clearData = buffer.remove(0, 1).remove(buffer.size() - 1, 1);

        switch (packageFeature)
        {
        case GetTestDataFromClient:
            emit readyMessageFromServer(QString("Data from: %1 : %2\n")
                                        .arg(camera->peerAddress().toString())
                                        .arg(QString(clearData)), camera);
            break;
        case GetBaseBallCoordinates:
            emit readyMessageFromServer(QString("Baseball coordinates!"), camera);
            break;
        case IsServerPrepareToGetStream:
            emit readyMessageFromServer(QString("Do we ready to get stream?"), camera);
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
        }

        }
    }
}


void CameraServer::appendFrameToBuffer(cv::Mat frame, QTcpSocket* camera)
{
    QMutexLocker lock(&streamMutex);
    cameras[camera].frameBuffer.append(frame);
    cameras[camera].lastFrame = frame;
    lock.unlock();
    emit frameFromStreamReady(camera);
}

void CameraServer::getFrameInternal(QTcpSocket* camera, qint32 port, const QString& path, bool sync = false)
{
    if (!cameras[camera].streamIsActive)
    {
        if (!cameras[camera].curParams.rawFrame)
        {
            //setenv("GST_DEBUG","2", 0);
            cv::VideoCapture cap(QString("udpsrc port=%1 ! application/x-rtp,media=video, payload=26,encoding-name=JPEG,framerate=30/1 !"
                                         " rtpjpegdepay ! jpegdec ! videoconvert ! appsink").arg(port).toStdString(),
                                 cv::CAP_GSTREAMER);
            cv::Mat frame;
            if (cap.isOpened())
            {
                while (true)
                {
                    cap.read(frame);
                    if (!frame.empty())
                    {
                        cv::Mat forSave;
                        cv::cvtColor(frame, forSave, CV_RGB2BGR);
                        if (!path.isEmpty())
                            qDebug() << cv::imwrite(path.toStdString(), forSave);
                        emit finishedGetData(camera);
                        appendFrameToBuffer(frame, camera);
                        break;
                    }
                }
            }
        }
        else
        {
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
                qDebug() << "new conn" << socket;
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
                                qDebug() << "INT TIME" << intTime;
                                emit syncTimeGot();
                                cameras[camera].syncTime = intTime;
                                cameras[camera].machineSyncTime = t;
                            }

                            qDebug()<< "read" << in.readRawData((char*)matRaw.data, width * height) << t << intTime;
                            buffer = buffer.remove(0, fullSize);
                            Mat mat;
                            cvtColor(matRaw, mat, CV_BayerBG2RGB);
                            appendFrameToBuffer(mat, camera);
                            if (!path.isEmpty())
                            {
                                Mat matBRG;
                                cvtColor(matRaw,matBRG, CV_BayerBG2BGR);
                                cv::imwrite(path.toStdString(), matBRG);
                            }
                            ++i;
                        }
                    }
                }
                server->close();
                //rawStreams.remove(server);
                //delete server;
            }
            else
            {
                qDebug() << "ERROR" << server->errorString();
            }
        }

    }
}

void CameraServer::writeRawVideoStream(QTcpSocket* camera, QByteArray& buffer, cv::Mat matRaw, cv::Mat mat, qint32& i, WriteRawVideoData &d)
{
    d.in.device()->reset();
    quint64 intTime, oldIntTime;
    d.in.readRawData((char*)&intTime, sizeof(quint64));
    oldIntTime = intTime;
    if (i == correctSyncTimeEvery)
    {
        cameras[camera].syncTime = intTime;
        correctSyncTime();
        i = 0;
    }
    if (cameras.firstKey() == camera)
    {
        intTime = intTime - syncTime;
    }
    qint32 exposureCenterOffset = (cameras[camera].curParams.exposure / 2) * exposureToCameraIntTime;
    intTime += exposureCenterOffset;
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
        cvtColor(matRaw, mat, CV_BayerBG2RGB);
        appendFrameToBuffer(mat, camera);
    }
    d.ints << intTime << "\t" << t.toString("hh:mm:ss.zzz") << "\t" << exposureCenterOffset << "\t" << syncTime << endl;
    ++i;
}

void CameraServer::getVideoInternal(QTcpSocket* camera, qint32 frameCount, qint32 port, bool compress, const QString& path = QString())
{
    try
    {
        if (!cameras[camera].streamIsActive)
        {
            if (!cameras[camera].curParams.rawFrame)
            {
                cameras[camera].streamIsActive = true;
                qDebug() << QString("udpsrc port=%1 ! application/x-rtp,encoding-name=H265,payload=96 ! rtph265depay "
                                    "! h265parse ! queue ! omxh265dec ! videoconvert ! appsink").arg(port);
                setenv("GST_DEBUG","2", 0);
                cv::VideoCapture cap(QString("udpsrc port=%1 ! application/x-rtp,encoding-name=H265,payload=96 ! rtph265depay "
                                             "! h265parse ! queue ! omxh265dec ! videoconvert ! appsink")
                                     .arg(port).toStdString(), cv::CAP_GSTREAMER);

                cv::Mat frame;
                if (cap.isOpened())
                {
                    qint32 i = 0;
                    while (true)
                    {
                        cap.read(frame);
                        if (!frame.empty())
                        {
                            if (cameras[camera].streamIsActive)
                            {
                                appendFrameToBuffer(frame, camera);
                            }
                            else
                            {
                                //сохранить видео
                            }
                            ++i;
                            if ((frameCount != -1 && i >= frameCount)
                                    || !cameras[camera].streamIsActive)
                            {
                                break;
                            }

                        }
                    }
                }
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
                    qDebug() << "new conn" << socket;
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
        data.append(EndOfMessageSign);
        socket->write(data);
    }
}

void CameraServer::requestCurrentCameraParameters(QTcpSocket* socket)
{
    QByteArray array;
    array.append(AskCurrentCameraParams);
    array.append(EndOfMessageSign);
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
        QtConcurrent::run(this, &CameraServer::getFrameInternal, socket, port, path, false);
        QByteArray data;
        data.append(RequestCameraFrame);
        data.append(EndOfMessageSign);
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
        data.append(EndOfMessageSign);
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
            data.append(EndOfMessageSign);
            socket->write(data);
            emit streamFromIsComing(socket);
            QtConcurrent::run(this, &CameraServer::getVideoInternal, socket, -1, port, false, QString());

        }
        else
        {
            data.append(StopStream);
            data.append(EndOfMessageSign);
            socket->write(data);
            cameras[socket].streamIsActive = false;

        }
    }
}

void CameraServer::restartCamera(QTcpSocket* socket)
{
    QByteArray data;
    data.append(RestartCamera);
    data.append(EndOfMessageSign);
    socket->write(data);
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

    for (auto& i : cameras)
    {
        delete i.timer;
    }
}
