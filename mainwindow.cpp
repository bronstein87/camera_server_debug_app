//#include <gst/gst.h>
//include <glib.h>
//#include <gst/rtsp-server/rtsp-server.h>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QImage>
#include <QPixmap>
#include <QVBoxLayout>
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
using namespace cv;


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), server(new CameraServer()),
    ui(new Ui::MainWindow), sync(new Synchronizer), optionsWindow(new CameraOptionsWindow(sync.data(), server)),
    serverOptionsWindow(new ServerOptionsWindow(sync.data()))
{

    ui->setupUi(this);
    connect(server, &CameraServer::readyMessageFromServer, this, [this](auto mes, auto socket)
    {
        if (!optionsWindow->showExposureVerbose() && mes.contains("AUTOEXP"))
        {
            return;
        }
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
        ui->cameraTabWidget->setCurrentIndex(ui->cameraTabWidget->count() - 1);
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
        Mat imgRgb;
        cvtColor(img, imgRgb, CV_BGR2RGB);
        static_cast <EditFrameQGraphicsScene*> (view->scene())
                ->setFrame(QPixmap::fromImage(
                               QImage(static_cast <uchar*>(imgRgb.data), img.cols, img.rows, img.step, QImage::Format_RGB888)));
    });


    connect(optionsWindow, &CameraOptionsWindow::correlatedImageReady, this, [this](auto mat, auto socket)
    {
        ZoomGraphicsView* view = camerasGUI[socket].view;
        Mat imgRgb;
        cvtColor(mat, imgRgb, CV_BGR2RGB);
        static_cast <EditFrameQGraphicsScene*> (view->scene())
                ->setFrame(QPixmap::fromImage(
                               QImage(static_cast <uchar*> (imgRgb.data), mat.cols, mat.rows, mat.step, QImage::Format_RGB888)));
    });



    connect(server, &CameraServer::autoCalibrateImageReady, this, [this](auto mat, auto socket)
    {
        if (optionsWindow->getShowAutoCalibrate())
        {
            ZoomGraphicsView* view = camerasGUI[socket].view;
            Mat imgRgb;
            cvtColor(mat, imgRgb, CV_BGR2RGB);
            static_cast <EditFrameQGraphicsScene*> (view->scene())
                    ->setFrame(QPixmap::fromImage(
                                   QImage(static_cast <uchar*> (imgRgb.data), mat.cols, mat.rows, mat.step, QImage::Format_RGB888)));
        }

    });

    connect(server, &CameraServer::resultPictureReady, this, [this](QImage& pic)
    {
        if (optionsWindow->getShowResults())
        {
            for (auto& i : camerasGUI.keys())
            {
                ZoomGraphicsView* view = camerasGUI[i].view;
                static_cast <EditFrameQGraphicsScene*> (view->scene())
                        ->setFrame(QPixmap::fromImage(pic));
            }
        }
    });

    QMap <qint32, QColor> plotMap;
    QVector <QCustomPlot*> plots;
    plotMap.insert(3850, QColor(255, 0, 0));
    plotMap.insert(4510, QColor(0, 255, 0));
    ApproximationVisualizer& av = ApproximationVisualizer::instance();
    av.initCamerasPlots(plotMap);
    av.setPlots(QVector <QCustomPlot*> {ui->errorsCustomPlot, ui->timesCustomPlot, ui->recCountCustomPlot, ui->calibDiffFirstCustomPlot,
                                        ui->calibDiffSecondCustomPlot});
    connect(&av, &ApproximationVisualizer::measureCountChanged, [this](qint32 count)
    {
        ui->plotHorizontalSlider->setRange(0, count);
    });

    QDir dir;
    QString gamesDirName = QString("games_%1").arg(QDate::currentDate().toString("dd_MM_yyyy"));
    if (dir.mkdir(gamesDirName))
    {
        qDebug() << dir.setCurrent(gamesDirName);
        dir.mkdir("4510");
        dir.mkdir("3850");
        dir.mkdir("results");
        dir.mkdir("graphs");
        dir.mkdir("pictures");
        dir.setCurrent(QApplication::applicationDirPath());
    }


    // av.drawTracerDebug(QString("D:/REC_CAMERAS/actual_server/build-camera_server_debug_app-Desktop_Qt_5_12_2_MSVC2017_64bit-Debug/games_17_09_2019/video3850_21_19_31.avi")
    //                   , QString("D:/REC_CAMERAS/actual_server/build-camera_server_debug_app-Desktop_Qt_5_12_2_MSVC2017_64bit-Debug/games_17_09_2019/video4510_21_19_31.avi"));
    //av.drawTracerDebug(QString("D:/REC_CAMERAS/actual_server/build-camera_server_debug_app-Desktop_Qt_5_12_2_MSVC2017_64bit-Debug/games_17_09_2019/video3850_20_56_13.avi")
    //                   , QString("D:/REC_CAMERAS/actual_server/build-camera_server_debug_app-Desktop_Qt_5_12_2_MSVC2017_64bit-Debug/games_17_09_2019/video4510_20_56_13.avi"));

    // rtsp://admin:dRX77QDV@10.20.55.110:554/cam/realmonitor?channel=1&subtype=0
    //    qputenv("GST_DEBUG", "4");
       // QString pipeLine = QString("rtspsrc location=rtsp://10.20.55.104:554/snl/live/1/1/Ux/sido=-QBP00YJW2yhJ ! rtph264depay ! h264parse ! avdec_h264 ! videoconvert ! appsink");
        QString pipeLine = QString("rtspsrc location=rtsp://admin:dRX77QDV@10.20.55.110:554/cam/realmonitor?channel=1&subtype=0 ! rtph264depay ! h264parse ! avdec_h264 ! videoconvert ! appsink");
        cv::VideoCapture cap(pipeLine.toStdString(), cv::CAP_GSTREAMER);
        if (cap.isOpened())
        {

            Mat m;
            double prev = 0;
            while (true)
            {
                if (cap.read(m))
                {
                    double cur = cap.get(CV_CAP_PROP_POS_MSEC);
                    qDebug() <<  m.cols << m.rows << cap.get(CV_CAP_PROP_POS_MSEC) << cur - prev;
                    prev = cur;
                }
            }
        }
    //    QElapsedTimer t;
    //    t.start();
    //    QThread::msleep(100);
    //    qDebug() << t.elapsed();
    //    QThread::msleep(100);
    //    qDebug() << t.elapsed();
//    QLinkedList <TestSt> ll;
//    for (qint32 i = 0; i < 10; ++i)
//    {
//        ll.append(TestSt());
//    }

//    for (qint32 j = 0; j < 10; ++j)
//    {
//        auto it = ll.begin();
//        while (it != ll.end())
//        {
//            --it->cnt;
//            ++it;
//        }
//    }
//    qDebug() << "qq";
//    BallApproximator approx;
//    QDir dirr("D:/REC_CAMERAS/actual_server/results/");
//    auto list = dirr.entryList(QDir::Files);
//    for (auto i : list)
//    {
//        server->testApproximation(QString("D:/REC_CAMERAS/actual_server/results/") + i, approx);
//    }

  // server->testApproximation("D:/REC_CAMERAS/actual_server/ideo_16_38_50", approx);
//    server->testApproximation("D:/REC_CAMERAS/actual_server/ideo_19_20_29", approx);
//    QImage img = av.makeShortPicture(approx, QString());
//    auto m = Mat(img.height(), img.width(), CV_8UC4, (void*)img.constBits());
//   cvtColor(m, m, CV_RGBA2BGR); // в отдельную функцию
//    HitParameters p;
//    p.angle = 10;
//    p.initSpeed = 80;
//    p.distance = 100;
//    auto img2 = av.addHitInfo(m, p);
//    img2.save("D:/REC_CAMERAS/test.jpg");

}

MainWindow::~MainWindow()
{
    delete server;
    delete optionsWindow;
    delete serverOptionsWindow;
    delete ui;
}




void MainWindow::on_startServerToolButton_clicked()
{
    server->startServer();
}

void MainWindow::on_cameraOptionsToolButton_clicked()
{
    if (ui->cameraTabWidget->count() >= 1)
    {
        optionsWindow->show();
    }

}

void MainWindow::on_serverOptionsToolButton_clicked()
{
    serverOptionsWindow->show();
}

void MainWindow::on_plotHorizontalSlider_sliderMoved(int position)
{
    ApproximationVisualizer& av = ApproximationVisualizer::instance();
    av.showMeasure(position, 3850, 4510); // tmp
}

void MainWindow::on_applyFirstCalibToolButton_clicked()
{
    CalibrationAdjustHelper helper;
    auto map = ApproximationVisualizer::instance().getCalibMap();
    helper.updateCalibrateInfo(3850, map[3850].last().eOr); //tmp
    imwrite(QString("calibrate/new%1.png").arg(3850).toStdString(), server->getLastCameraFrame(3850));
}

void MainWindow::on_applySecondCalibToolButton_clicked()
{
    CalibrationAdjustHelper helper;
    auto map = ApproximationVisualizer::instance().getCalibMap();
    helper.updateCalibrateInfo(4510, map[4510].last().eOr); //tmp
    imwrite(QString("calibrate/new%1.png").arg(4510).toStdString(), server->getLastCameraFrame(4510));
}
