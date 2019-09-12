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
#include <approximationvisualizer.h>

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


private slots:


    void on_startServerToolButton_clicked();

    void on_cameraOptionsToolButton_clicked();

    void on_serverOptionsToolButton_clicked();

    void on_plotHorizontalSlider_sliderMoved(int position);

    void on_applyFirstCalibToolButton_clicked();

    void on_applySecondCalibToolButton_clicked();

private:

    QMap <QTcpSocket*, CameraGUIInfo> camerasGUI;
    CameraServer* server;
    Ui::MainWindow *ui;
    QScopedPointer <Synchronizer> sync;
    CameraOptionsWindow* optionsWindow = nullptr;
    ServerOptionsWindow* serverOptionsWindow = nullptr;
    cv::Mat img;
    QImage* imgIn;
    CalibrationAdjustHelper helper;




};

#endif // MAINWINDOW_H
