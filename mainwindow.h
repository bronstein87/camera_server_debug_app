#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <cameraserver.h>
#include <editframeqgraphicsview.h>
#include <cameraserver.h>
#include <zoomgraphicsview.h>
#include <QPlainTextEdit>
#include <QLabel>
#include <QMap>
#include <cameraoptionswindow.h>
#include <serveroptionswindow.h>
#include <calibrationadjusthelper.h>

class CameraOptionsWindow;
struct CameraGUIInfo
{
    ZoomGraphicsView* view;
    QPlainTextEdit* infoWindow;
    qint32 tabIndex;
};

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

signals:
    //currentCameraChanged(QTcpSocket* cam);

private slots:


    void on_serverStartAction_triggered();

    void on_openOptionsAction_triggered();

    void on_serverOptionsAction_triggered();


    void on_pushButton_clicked();

private:

    QMap <QTcpSocket*, CameraGUIInfo> camerasGUI;
    CameraServer* server;
    Ui::MainWindow *ui;
    QScopedPointer <Synchronizer> sync;
    CameraOptionsWindow* optionsWindow = nullptr;
    ServerOptionsWindow* serverOptionsWindow = nullptr;
    EditFrameQGraphicsScene scene;
    cv::Mat img;
    QImage* imgIn;
    CalibrationAdjustHelper helper;



};

#endif // MAINWINDOW_H
