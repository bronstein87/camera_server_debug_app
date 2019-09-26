#ifndef CAMERAOPTIONSWINDOW_H
#define CAMERAOPTIONSWINDOW_H

#include <QWidget>
#include <cameraserver.h>
#include <mainwindow.h>
#include <QFileDialog>
#include <QMessageBox>
#include <synchronizer.h>
#include <calibrationadjusthelper.h>
#include <QSettings>


class MainWindow;

namespace Ui {
class CameraOptionsWindow;
}

class CameraOptionsWindow : public QWidget
{
    Q_OBJECT

public:
    explicit CameraOptionsWindow(Synchronizer* snc, CameraServer* server, QWidget* parent = 0);
    ~CameraOptionsWindow();

    void setCurrentCamera(QTcpSocket* socket, EditFrameQGraphicsScene* scene);

    bool showExposureVerbose();

    bool getShowResults();

    bool getShowAutoCalibrate();

private slots:


    void on_focusingSpinBox_editingFinished();

    void on_pictureParamGroupBox_toggled(bool arg1);

    void on_requestFramePushButton_clicked();

    void on_requestVideoPushButton_clicked();

    void on_gammaSlider_valueChanged(int value);

    void on_gainSlider_valueChanged(int value);

    void on_sharpSlider_valueChanged(int value);

    void on_sendFrameRateSpinBox_editingFinished();

    void on_startStreamPushButton_clicked();

    void on_stopStreamPushButton_clicked();

    void on_chooseFramePathToolButton_clicked();

    void on_chooseVideoPathToolButton_clicked();

    void on_gammaLineEdit_editingFinished();

    void on_gainLineEdit_editingFinished();

    void on_sharpLineEdit_editingFinished();

    void on_whiteBalanceCheckBox_clicked(bool checked);

    void on_equalizationCheckBox_clicked(bool checked);

    void on_streamPortSpinBox_editingFinished();

    void on_exposureHorizontalSlider_valueChanged(int value);

    void on_pixelClockHorizontalSlider_valueChanged(int value);

    void on_frameRateHorizontalSlider_valueChanged(int value);

    void on_autoMinGainSpinBox_editingFinished();

    void on_autoPercentDoubleSpinBox_editingFinished();

    void on_enableAutoExposureCheckBox_clicked(bool checked);

    void on_dayRadioButton_clicked(bool checked);

    void on_eveningRadioButton_clicked(bool checked);

    void on_saturationHorizontalSlider_valueChanged(int value);

    void on_saturationLineEdit_editingFinished();

    void on_saturationRHorizontalSlider_valueChanged(int value);

    void on_saturationRLineEdit_editingFinished();

    void on_saturationGHorizontalSlider_valueChanged(int value);

    void on_saturationBLineEdit_editingFinished();

    void on_saturationBHorizontalSlider_valueChanged(int value);

    void on_shadowCoefDoubleSpinBox_editingFinished();

    void on_shadowGaussWindowSizeSpinBox_editingFinished();

    void on_hueHorizontalSlider_valueChanged(int value);

    void on_hueLineEdit_editingFinished();

    void on_autoLeftBorderDoubleSpinBox_editingFinished();

    void on_autoRightBorderDoubleSpinBox_editingFinished();

    void on_saturationGLineEdit_editingFinished();

    void on_autoMaxGainSpinBox_editingFinished();

    void on_shadowThresholdDoubleSpinBox_editingFinished();

    void on_chooseWBROIToolButton_clicked();

    void on_showWBROIToolButton_clicked();

    void on_chooseMainROIPushButton_clicked();

    void on_showMainROIPushButton_clicked();

    void on_chooseFirstROIPushButton_clicked();

    void on_chooseSecondROIPushButton_clicked();

    void on_showSecondROIPushButton_clicked();

    void on_showFirstROIPushButton_clicked();

    void on_rawTcpStreamCheckBox_clicked(bool checked);

    void on_durationSpinBox_editingFinished();

    void on_triggerModeGroupBox_clicked(bool checked);

    void on_highLowRadioButton_clicked(bool checked);

    void on_lowHighRadioButton_clicked(bool checked);

    void on_turnOnAutoVideoSaveRadioButton_clicked(bool checked);

    void on_turnOffAutoVideoSaveRadioButton_clicked(bool checked);

    void on_chooseDirectoryToolButton_clicked();

    void on_corrDoubleSpinBox_editingFinished();

    void on_searchAreaSizeSpinBox_editingFinished();

    void on_minSkoTemplateDoubleSpinBox_editingFinished();

    void on_maxAngleBetwDirDoubleSpinBox_editingFinished();

    void on_minSpeedDoubleSpinBox_editingFinished();

    void on_cannyMinSpinBox_editingFinished();

    void on_cannyMaxSpinBox_editingFinished();

    void on_minAreaSpinBox_editingFinished();

    void on_maxAreaSpinBox_editingFinished();

    void on_circularityDoubleSpinBox_editingFinished();

    void on_enableRecognizeCheckBox_toggled(bool checked);

    void on_ballRecognizeStepSpinBox_editingFinished();

    void on_debugRecFlagCheckBox_clicked(bool checked);

    void on_checkSynchroPushButton_clicked();

    void on_checkSyncPushButton_clicked();

    void on_restartCameraPushButton_clicked();

    void on_showVideoCheckBox_clicked(bool checked);

    void on_saveVideoCheckBox_clicked(bool checked);

    void on_syncFramePushButton_clicked();

    void on_chooseReferenceFramePushButton_clicked();

    void on_showPushButton_clicked();

    void on_calibratePushButton_clicked();

    void on_debounceEnableGroupBox_toggled(bool arg1);

    void on_debounceSpinBox_editingFinished();

    void on_saveCameraSettingsPushButton_clicked();

    void on_skoCountDoubleSpinBox_editingFinished();


    void on_fullPictureCheckBox_clicked(bool checked);

    void on_shortPictureCheckBox_clicked(bool checked);

    void on_turnOnSyncFrameEveryCheckBox_toggled(bool checked);

    void on_updateGraphsCheckBox_toggled(bool checked);

    void on_turnOnRecogAllCamerasCheckBox_toggled(bool checked);

    void on_compareCurrentRadioButton_toggled(bool checked);

    void on_compareReferenceRadioButton_toggled(bool checked);

    void on_compareNeibRadioButton_toggled(bool checked);

    void on_xDirDoubleSpinBox_editingFinished();

    void on_yDirDoubleSpinBox_editingFinished();

    void on_turnOnScenarioCheckBox_toggled(bool checked);

    void on_chooseMainHitROIPushButton_clicked();

    void on_showMainHitPushButton_clicked();

signals:
    void correlatedImageReady(cv::Mat img, QTcpSocket* socket);

private:

    void loadSettings();

    void saveSettings();

    void handleCurrentCameraParams(QTcpSocket *camera, const CurrentCameraParams& params);

    Ui::CameraOptionsWindow *ui;
    CameraServer* cameraServer;
    Synchronizer* sync = nullptr;
    QTcpSocket* currentCamera = nullptr;
    EditFrameQGraphicsScene* roiScene;
    CalibrationAdjustHelper calibAdjustHelper;
    CompareFlag flag = CompareFlag::None;



};

#endif // CAMERAOPTIONSWINDOW_H
