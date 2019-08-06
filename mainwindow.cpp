#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QImage>
#include <QPixmap>
#include <QVBoxLayout>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), server(new CameraServer(this)),
    ui(new Ui::MainWindow), sync(new Synchronizer), optionsWindow(new CameraOptionsWindow(sync.data(), server)),  serverOptionsWindow(new ServerOptionsWindow(sync.data()))
{

    ui->setupUi(this);
    connect(server, &CameraServer::readyMessageFromServer, this, [this](auto mes, auto socket)
    {
        QString currentDt = QDateTime::currentDateTime().toString(Qt::ISODate);
        if (socket != nullptr)
        {
            camerasGUI[socket].infoWindow->appendPlainText(QString("%1 %2")
                                                           .arg(currentDt)
                                                           .arg(mes));
        }
        else
        {
            qDebug() << mes;
        }

    });


    connect(ui->cameraTabWidget, &QTabWidget::currentChanged, [this](auto index)
    {
        QMapIterator <QTcpSocket*, CameraGUIInfo> i(camerasGUI);
        while (i.hasNext())
        {
            i.next();
            if (i.value().tabIndex == index)
            {
                optionsWindow->setCurrentCamera(i.key(), static_cast <EditFrameQGraphicsScene*> (i.value().view->scene()));
                server->requestCurrentCameraParameters(i.key());
            }
        }
    });

    connect(server, &CameraServer::connectedCameraSocket, this, [this] (QTcpSocket* socket)
    {
        QWidget* wgt = new QWidget(ui->cameraTabWidget);
        QVBoxLayout* layout = new QVBoxLayout(wgt);
        auto view = new ZoomGraphicsView(wgt);
        auto monitor = new QPlainTextEdit(wgt);
        monitor->setMaximumHeight(200);
        CameraGUIInfo info;
        info.infoWindow = monitor;
        info.view = view;
        view->setScene(new EditFrameQGraphicsScene(view));

        layout->addWidget(view);
        layout->addWidget(new QLabel("Мониторинг", wgt));
        layout->addWidget(monitor);
        wgt->setLayout(layout);
        info.tabIndex = ui->cameraTabWidget->addTab(wgt, QString("Камера №%1").arg(camerasGUI.size() + 1));

        camerasGUI.insert(socket, info);
        optionsWindow->setCurrentCamera(socket, static_cast <EditFrameQGraphicsScene*> (view->scene()));
        server->requestCurrentCameraParameters(socket);
    });

    connect(server, &CameraServer::disconnectedCameraSocket, this, [this] (QTcpSocket* socket)
    {
        qint32 index = camerasGUI[socket].tabIndex;
        camerasGUI.remove(socket);
        ui->cameraTabWidget->removeTab(index);
        if (ui->cameraTabWidget->count() == 0)
        {
            optionsWindow->hide();
        }
    });

    connect(server, &CameraServer::frameFromStreamReady, this, [this](auto socket)
    {
        ZoomGraphicsView* view = camerasGUI[socket].view;
        auto img = server->getFrameFromStream(socket);
        static_cast <EditFrameQGraphicsScene*> (view->scene())
                ->setFrame(QPixmap::fromImage(
                               QImage((uchar*)img.data, img.cols, img.rows, img.step, QImage::Format_RGB888)));
    });


    connect(optionsWindow, &CameraOptionsWindow::correlatedImageReady, this, [this](auto mat, auto socket)
    {
        ZoomGraphicsView* view = camerasGUI[socket].view;
        static_cast <EditFrameQGraphicsScene*> (view->scene())
                ->setFrame(QPixmap::fromImage(
                               QImage((uchar*)mat.data, mat.cols, mat.rows, mat.step, QImage::Format_RGB888)));
    });


    QPixmap px("/home/nvidia/actual_server/server_debug/calibrate/picture.jpg");
    QGraphicsPixmapItem* item  = new QGraphicsPixmapItem(px);

     server->testApproximation("/home/nvidia/actual_server/server_debug/calibrate/4510_approx",
                               "/home/nvidia/actual_server/server_debug/calibrate/3850_approx", item);
     QGraphicsScene* scene = new QGraphicsScene(this);
     scene->addItem(item);
     ui->graphicsView->setScene(scene);
     scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
     QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
     image.fill(Qt::transparent);                                              // Start all pixels transparent

     QPainter painter(&image);
     scene->render(&painter);
     image.save("/home/nvidia/actual_server/server_debug/calibrate/file_name.png");






   //server->createCalibrateImage("/home/nvidia/recalibrate_debug/right (3850).bmp.tif", "/home/nvidia/recalibrate_debug/3850_wh", "/home/nvidia/recalibrate_debug/cor2.png", 35);
   // Mat fm(1216, 1936, CV_8UC3);
   // Mat m = imread("/home/nvidia/actual_server/server_debug/calibrate/4510test.bmp");
    //m.copyTo(fm(Rect(0, 0, 1920, 1080)));
   // imwrite("/home/nvidia/actual_server/server_debug/calibrate/TEST.png", fm);
   // helper.createCalibrateImage("/home/nvidia/actual_server/server_debug/calibrate/left (4510).bmp.tif",
   //                           "/home/nvidia/actual_server/server_debug/calibrate/4510_wh", fm, 25);

//    Calibration::ExteriorOr eOr, nEOr;
//    Calibration::SpacecraftPlatform::CAMERA::CameraParams cam, nCam;
//    CalibrationAdjustHelper::readCurrentCalibrationParameters(4510, "/home/nvidia/actual_server/server_debug/calibrate", eOr, cam);
//    helper.recalibrate(QVector <qint32> {}, eOr, cam, nEOr, nCam);
//    qDebug() << "finished" << nEOr.Point.X << nEOr.Point.Y << nEOr.Point.Z
//             << Calibration::Radian2Deg(nEOr.OPK.Omega)
//             << Calibration::Radian2Deg(nEOr.OPK.Phi)
//             << Calibration::Radian2Deg(nEOr.OPK.Kappa);
}

MainWindow::~MainWindow()
{
    delete ui;
}



//    imgIn = new QImage((uchar*) img.data, img.cols, img.rows, img.step, QImage::Format_RGB888);

void MainWindow::on_serverStartAction_triggered()
{
    server->startServer();
}

void MainWindow::on_openOptionsAction_triggered()
{
    if (ui->cameraTabWidget->count() >= 1)
    {
        optionsWindow->show();
    }
}

void MainWindow::on_serverOptionsAction_triggered()
{
    serverOptionsWindow->show();
}
#include <unistd.h>




void MainWindow::on_pushButton_clicked()
{
//    sync->blockSignals(true);
//    if (sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst, 1))
//    {
//        usleep(1000000);
//        if (sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x1))
//        {
//            //cameraServer->checkSync();
//            QTimer::singleShot(750, this, [this]()
//            {
//                sync->sendCommand(Synchronizer::SynchronizerProtocol::StopSync, 0x1);
//                usleep(1000000);
//            }
//            );

//        }
//    }
//    sync->blockSignals(false);
    sync->blockSignals(true);
    if (sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst, 1))
    {
        if (sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x1))
        {
            //cameraServer->checkSync();
            QTimer::singleShot(500, this, [this]()
            {
                sync->sendCommand(Synchronizer::SynchronizerProtocol::StopSync, 0x1);
            }
            );

        }
    }
    sync->blockSignals(false);
}
