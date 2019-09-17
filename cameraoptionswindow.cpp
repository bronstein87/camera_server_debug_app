#include "cameraoptionswindow.h"
#include "ui_cameraoptionswindow.h"
//#include <unistd.h>
using namespace cv;

#define IS_SET_TRIGGER_OFF                  0x0000
#define IS_SET_TRIGGER_CONTINUOUS           0x1000
#define IS_SET_TRIGGER_HI_LO                (IS_SET_TRIGGER_CONTINUOUS | 0x0001)
#define IS_SET_TRIGGER_LO_HI                (IS_SET_TRIGGER_CONTINUOUS | 0x0002)

CameraOptionsWindow::CameraOptionsWindow(Synchronizer* snc, CameraServer* server, QWidget* parent) :
    cameraServer(server), sync(snc), QWidget(parent),
    ui(new Ui::CameraOptionsWindow)
{
    ui->setupUi(this);
    setWindowFlags(Qt::WindowStaysOnTopHint);
    connect(cameraServer, &CameraServer::currentCameraParamsReady, this, &CameraOptionsWindow::handleCurrentCameraParams);
    ui->rawTcpStreamCheckBox->setChecked(false);
    connect(sync, &Synchronizer::changed, cameraServer, &CameraServer::updateParamsForAllCameras);
    loadSettings();
}



void CameraOptionsWindow::handleCurrentCameraParams(QTcpSocket* camera, const CurrentCameraParams& params)
{
    if (camera == currentCamera)
    {
        ui->exposureHorizontalSlider->blockSignals(true);

        ui->exposureHorizontalSlider->setMinimum(params.minExposure * 100);
        ui->minExpLabel->setText(QString::number(params.minExposure));

        ui->exposureHorizontalSlider->setMaximum(params.maxExposure * 100);
        ui->maxExpLabel->setText(QString::number(params.maxExposure));

        ui->exposureHorizontalSlider->setValue(params.exposure * 100);
        ui->curExpLabel->setText(QString::number(params.exposure));
        ui->exposureHorizontalSlider->blockSignals(false);

        ui->pixelClockHorizontalSlider->blockSignals(true);

        ui->pixelClockHorizontalSlider->setMinimum(params.minPixelClock);
        ui->minPixelClockLabel->setText(QString::number(ui->pixelClockHorizontalSlider->minimum()));

        ui->pixelClockHorizontalSlider->setMaximum(params.maxPixelClock);
        ui->maxPixelClockLabel->setText(QString::number(ui->pixelClockHorizontalSlider->maximum()));

        ui->pixelClockHorizontalSlider->setValue(params.pixelClock);
        ui->curPixelClockLabel->setText(QString::number(ui->pixelClockHorizontalSlider->value()));
        ui->pixelClockHorizontalSlider->blockSignals(false);


        ui->frameRateHorizontalSlider->blockSignals(true);

        ui->frameRateHorizontalSlider->setMinimum(params.minFrameRate);
        ui->minFrameRateLabel->setText(QString::number(ui->frameRateHorizontalSlider->minimum()));

        ui->frameRateHorizontalSlider->setMaximum(params.maxFrameRate);
        ui->maxFrameRateLabel->setText(QString::number(ui->frameRateHorizontalSlider->maximum()));

        ui->frameRateHorizontalSlider->setValue(params.frameRate);
        ui->curFrameRateLabel->setText(QString::number(ui->frameRateHorizontalSlider->value()));
        ui->frameRateHorizontalSlider->blockSignals(false);

        ui->sendFrameRateSpinBox->blockSignals(true);
        ui->sendFrameRateSpinBox->setValue(params.sendFrameRate);
        ui->sendFrameRateSpinBox->blockSignals(false);
        ui->dayRadioButton->blockSignals(true);
        ui->dayRadioButton->blockSignals(true);
        if (params.autoExpParams.light == Day)
        {
            ui->dayRadioButton->setChecked(true);
        }
        else
        {
            ui->eveningRadioButton->setChecked(true);
        }
        ui->dayRadioButton->blockSignals(false);
        ui->dayRadioButton->blockSignals(false);

        ui->autoMinGainSpinBox->blockSignals(true);
        ui->autoMinGainSpinBox->setValue(params.autoExpParams.minGainCoeff);
        ui->autoMinGainSpinBox->blockSignals(false);

        ui->autoMaxGainSpinBox->blockSignals(true);
        ui->autoMaxGainSpinBox->setValue(params.autoExpParams.maxGainCoeff);
        ui->autoMaxGainSpinBox->blockSignals(false);

        ui->autoPercentDoubleSpinBox->blockSignals(true);
        ui->autoPercentDoubleSpinBox->setValue(params.autoExpParams.maxPercent * 100);
        ui->autoPercentDoubleSpinBox->blockSignals(false);

        ui->autoLeftBorderDoubleSpinBox->blockSignals(true);
        ui->autoLeftBorderDoubleSpinBox->setValue(params.autoExpParams.minRelCoef);
        ui->autoLeftBorderDoubleSpinBox->blockSignals(false);

        ui->autoRightBorderDoubleSpinBox->blockSignals(true);
        ui->autoRightBorderDoubleSpinBox->setValue(params.autoExpParams.maxRelCoef);
        ui->autoRightBorderDoubleSpinBox->blockSignals(false);

        ui->enableAutoExposureCheckBox->blockSignals(true);
        ui->enableAutoExposureCheckBox->setChecked(params.autoExposureIntenal);
        ui->enableAutoExposureCheckBox->blockSignals(false);

        ui->gainSlider->blockSignals(true);
        ui->gainSlider->setValue(params.gain);
        ui->gainSlider->blockSignals(false);

        ui->gainLineEdit->blockSignals(true);
        ui->gainLineEdit->setText(QString::number(params.gain));
        ui->gainLineEdit->blockSignals(false);


        ui->rawTcpStreamCheckBox->blockSignals(true);
        ui->rawTcpStreamCheckBox->setChecked(params.rawFrame);
        ui->rawTcpStreamCheckBox->blockSignals(false);

        ui->durationSpinBox->blockSignals(true);
        ui->durationSpinBox->setValue(params.videoDuration);
        ui->durationSpinBox->blockSignals(false);

        ui->widthSpinBox->blockSignals(true);
        ui->widthSpinBox->setValue(params.width);
        ui->widthSpinBox->blockSignals(false);

        ui->heightSpinBox->blockSignals(true);
        ui->heightSpinBox->setValue(params.height);
        ui->heightSpinBox->blockSignals(false);

        ui->triggerModeGroupBox->blockSignals(true);
        ui->triggerModeGroupBox->setChecked(params.triggerMode);
        ui->highLowRadioButton->blockSignals(true);
        ui->lowHighRadioButton->blockSignals(true);
        if (params.triggerMode == IS_SET_TRIGGER_HI_LO)
        {
            ui->highLowRadioButton->setChecked(true);
        }
        else if (params.triggerMode == IS_SET_TRIGGER_LO_HI)
        {
            ui->lowHighRadioButton->setChecked(true);
        }
        ui->highLowRadioButton->blockSignals(false);
        ui->lowHighRadioButton->blockSignals(false);
        ui->triggerModeGroupBox->blockSignals(false);

        ui->pictureParamGroupBox->blockSignals(true);
        ui->pictureParamGroupBox->setChecked(params.pictureParamFlag);
        ui->pictureParamGroupBox->blockSignals(false);

        ui->whiteBalanceCheckBox->blockSignals(true);
        ui->whiteBalanceCheckBox->setChecked(params.whiteBalance);
        ui->whiteBalanceCheckBox->blockSignals(false);

        ui->equalizationCheckBox->blockSignals(true);
        ui->equalizationCheckBox->setChecked(params.equalization);
        ui->equalizationCheckBox->blockSignals(false);

        ui->enableAutoExposureCheckBox->blockSignals(true);
        ui->enableAutoExposureCheckBox->setChecked(params.autoExposureIntenal);
        ui->enableAutoExposureCheckBox->blockSignals(false);

        ui->focusingSpinBox->blockSignals(true);
        ui->focusingSpinBox->setValue(params.focusing);
        ui->focusingSpinBox->blockSignals(false);

        ui->sendFrameRateSpinBox->blockSignals(true);
        ui->sendFrameRateSpinBox->setValue(params.sendFrameRate);
        ui->sendFrameRateSpinBox->blockSignals(false);

        ui->gammaSlider->blockSignals(true);
        ui->gammaLineEdit->blockSignals(true);
        ui->gammaSlider->setValue(params.gamma * 100);
        ui->gammaLineEdit->setText(QString::number(params.gamma));
        ui->gammaSlider->blockSignals(false);
        ui->gammaLineEdit->blockSignals(false);

        ui->sharpSlider->blockSignals(true);
        ui->sharpLineEdit->blockSignals(true);
        ui->sharpSlider->setValue(params.sharp * 10);
        ui->sharpLineEdit->setText(QString::number(params.sharp));
        ui->sharpSlider->blockSignals(false);
        ui->sharpLineEdit->blockSignals(false);

        ui->saturationHorizontalSlider->blockSignals(true);
        ui->saturationLineEdit->blockSignals(true);
        ui->saturationHorizontalSlider->setValue(params.saturation);
        ui->saturationLineEdit->setText(QString::number(params.saturation));
        ui->saturationHorizontalSlider->blockSignals(false);
        ui->saturationLineEdit->blockSignals(false);


        ui->saturationRHorizontalSlider->blockSignals(true);
        ui->saturationRLineEdit->blockSignals(true);
        ui->saturationRHorizontalSlider->setValue(params.rSaturation);
        ui->saturationRLineEdit->setText(QString::number(params.rSaturation));
        ui->saturationRHorizontalSlider->blockSignals(false);
        ui->saturationRLineEdit->blockSignals(false);


        ui->saturationGHorizontalSlider->blockSignals(true);
        ui->saturationGLineEdit->blockSignals(true);
        ui->saturationGHorizontalSlider->setValue(params.gSaturation);
        ui->saturationGLineEdit->setText(QString::number(params.gSaturation));
        ui->saturationGHorizontalSlider->blockSignals(false);
        ui->saturationGLineEdit->blockSignals(false);


        ui->saturationBHorizontalSlider->blockSignals(true);
        ui->saturationBLineEdit->blockSignals(true);
        ui->saturationBHorizontalSlider->setValue(params.bSaturation);
        ui->saturationBLineEdit->setText(QString::number(params.bSaturation));
        ui->saturationBHorizontalSlider->blockSignals(false);
        ui->saturationBLineEdit->blockSignals(false);

        ui->hueHorizontalSlider->blockSignals(true);
        ui->hueLineEdit->blockSignals(true);
        ui->hueHorizontalSlider->setValue(params.hue);
        ui->hueLineEdit->setText(QString::number(params.hue));
        ui->hueHorizontalSlider->blockSignals(false);
        ui->hueLineEdit->blockSignals(false);

        ui->shadowThresholdDoubleSpinBox->blockSignals(true);
        ui->shadowThresholdDoubleSpinBox->setValue(params.shadowThreshold);
        ui->shadowThresholdDoubleSpinBox->blockSignals(false);

        ui->shadowCoefDoubleSpinBox->blockSignals(true);
        ui->shadowCoefDoubleSpinBox->setValue(params.shadowCoef);
        ui->shadowCoefDoubleSpinBox->blockSignals(false);

        ui->streamPortSpinBox->blockSignals(true);
        ui->streamPortSpinBox->setValue(params.portSendStream);
        ui->streamPortSpinBox->blockSignals(false);

        ui->shadowGaussWindowSizeSpinBox->blockSignals(true);
        ui->shadowGaussWindowSizeSpinBox->setValue(params.shadowGaussWindowSize);
        ui->shadowGaussWindowSizeSpinBox->blockSignals(false);

        ui->enableRecognizeCheckBox->blockSignals(true);
        ui->enableRecognizeCheckBox->setChecked(params.ballRecognizeFlag);
        ui->enableRecognizeCheckBox->blockSignals(false);

        ui->ballRecognizeStepSpinBox->blockSignals(true);
        ui->ballRecognizeStepSpinBox->setValue(params.ballRecognizeStep);
        ui->ballRecognizeStepSpinBox->blockSignals(false);

        ui->corrDoubleSpinBox->blockSignals(true);
        ui->corrDoubleSpinBox->setValue(params.recParams.corrCoef);
        ui->corrDoubleSpinBox->blockSignals(false);

        ui->skoCountDoubleSpinBox->blockSignals(true);
        ui->skoCountDoubleSpinBox->setValue(params.recParams.skoCoef);
        ui->skoCountDoubleSpinBox->blockSignals(false);

        ui->searchAreaSizeSpinBox->blockSignals(true);
        ui->searchAreaSizeSpinBox->setValue(params.recParams.searchAreaSize);
        ui->searchAreaSizeSpinBox->blockSignals(false);

        ui->minAreaSpinBox->blockSignals(true);
        ui->minAreaSpinBox->setValue(params.recParams.minArea);
        ui->minAreaSpinBox->blockSignals(false);

        ui->maxAreaSpinBox->blockSignals(true);
        ui->maxAreaSpinBox->setValue(params.recParams.maxArea);
        ui->maxAreaSpinBox->blockSignals(false);

        ui->minSkoTemplateDoubleSpinBox->blockSignals(true);
        ui->minSkoTemplateDoubleSpinBox->setValue(params.recParams.minSkoOnTemplate);
        ui->minSkoTemplateDoubleSpinBox->blockSignals(false);

        ui->maxAngleBetwDirDoubleSpinBox->blockSignals(true);
        ui->maxAngleBetwDirDoubleSpinBox->setValue(params.recParams.maxAngleBetwDirections);
        ui->maxAngleBetwDirDoubleSpinBox->blockSignals(false);

        ui->minSpeedDoubleSpinBox->blockSignals(true);
        ui->minSpeedDoubleSpinBox->setValue(params.recParams.minSpeed);
        ui->minSpeedDoubleSpinBox->blockSignals(false);

        ui->cannyMinSpinBox->blockSignals(true);
        ui->cannyMinSpinBox->setValue(params.recParams.cannyThresMin);
        ui->cannyMinSpinBox->blockSignals(false);

        ui->cannyMaxSpinBox->blockSignals(true);
        ui->cannyMaxSpinBox->setValue(params.recParams.cannyThresMax);
        ui->cannyMaxSpinBox->blockSignals(false);

        ui->circularityDoubleSpinBox->blockSignals(true);
        ui->circularityDoubleSpinBox->setValue(params.recParams.circularityCoeff);
        ui->circularityDoubleSpinBox->blockSignals(false);
        //RecognizeParameters recParams;

        ui->debugRecFlagCheckBox->blockSignals(true);
        ui->debugRecFlagCheckBox->setChecked(params.debugRecFlag);
        ui->debugRecFlagCheckBox->blockSignals(false);

        ui->showVideoCheckBox->blockSignals(true);
        ui->showVideoCheckBox->setChecked(cameraServer->getShowVideo(currentCamera));
        ui->showVideoCheckBox->blockSignals(false);

        ui->saveVideoCheckBox->blockSignals(true);
        ui->saveVideoCheckBox->setChecked(cameraServer->getWriteVideo(currentCamera));
        ui->saveVideoCheckBox->blockSignals(false);
        ui->debounceEnableGroupBox->blockSignals(true);
        ui->debounceEnableGroupBox->setChecked(params.debounceEnable);
        ui->debounceEnableGroupBox->blockSignals(false);
        ui->debounceSpinBox->blockSignals(true);
        ui->debounceSpinBox->setValue(params.debounceValue);
        ui->debounceSpinBox->blockSignals(false);
    }

}

CameraOptionsWindow::~CameraOptionsWindow()
{
    saveSettings();
    delete ui;
}

void CameraOptionsWindow::setCurrentCamera(QTcpSocket* socket, EditFrameQGraphicsScene* scene)
{
    currentCamera = socket;
    roiScene = scene;
    setWindowTitle(QString("Настройки камеры (%1)").arg(cameraServer->getLastCurrentParameters(socket).portSendStream));
}

bool CameraOptionsWindow::showExposureVerbose()
{
    return ui->exposureMessageCheckBox->isChecked();
}

void CameraOptionsWindow::on_gammaSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.gamma = (double)value / 100;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->gammaLineEdit->setText(QString::number(opt.gamma));
}

void CameraOptionsWindow::on_gainSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.gain = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->gainLineEdit->setText(QString::number(value));
}

void CameraOptionsWindow::on_sharpSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.sharp = (double)value / 10;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->sharpLineEdit->setText(QString::number(opt.sharp));
}




void CameraOptionsWindow::on_focusingSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.focusing = ui->focusingSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}


void CameraOptionsWindow::on_pictureParamGroupBox_toggled(bool arg1)
{
    CameraOptions opt;
    opt.pictureParamFlag = arg1;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_requestFramePushButton_clicked()
{
    //    if (ui->framePathLineEdit->text().isEmpty())
    //    {
    //        QMessageBox::warning(this, "Внимание", "Не указан путь!");
    //        return;
    //    }
    cameraServer->requestFrameFromCamera(currentCamera, ui->streamPortSpinBox->value(), ui->framePathLineEdit->text());
}

void CameraOptionsWindow::on_requestVideoPushButton_clicked()
{
    //    if (ui->videoPathLineEdit->text().isEmpty())
    //    {
    //        QMessageBox::warning(this, "Внимание", "Не указан путь!");
    //        return;
    //    }
    cameraServer->requestVideoFromCamera(ui->streamPortSpinBox->value(), ui->durationSpinBox->value(), ui->compressCheckBox->isChecked(), currentCamera,  ui->videoPathLineEdit->text());
}



void CameraOptionsWindow::on_sendFrameRateSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.frameRateSend = ui->sendFrameRateSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_startStreamPushButton_clicked()
{
    cameraServer->startStream(ui->streamPortSpinBox->value(), currentCamera, true);
}

void CameraOptionsWindow::on_stopStreamPushButton_clicked()
{
    cameraServer->startStream(-1, currentCamera, false);
}

void CameraOptionsWindow::on_chooseFramePathToolButton_clicked()
{
    ui->framePathLineEdit->setText(QFileDialog::getSaveFileName(this,
                                                                tr("Save Image"), ui->framePathLineEdit->text(), tr("Image Files (*.bmp)")));
}

void CameraOptionsWindow::on_chooseVideoPathToolButton_clicked()
{
    ui->videoPathLineEdit->setText(QFileDialog::getSaveFileName(this,
                                                                tr("Save video"), ui->videoPathLineEdit->text(), tr("Image Files (*.avi)")));
}



void CameraOptionsWindow::on_gammaLineEdit_editingFinished()
{
    ui->gammaSlider->setValue(ui->gammaLineEdit->text().toDouble() * 10);
}

void CameraOptionsWindow::on_gainLineEdit_editingFinished()
{
    ui->gainSlider->setValue(ui->gainLineEdit->text().toDouble());
}

void CameraOptionsWindow::on_sharpLineEdit_editingFinished()
{
    ui->sharpSlider->setValue(ui->sharpLineEdit->text().toDouble() * 10);
}

void CameraOptionsWindow::on_whiteBalanceCheckBox_clicked(bool checked)
{
    CameraOptions opt;
    opt.whiteBalance = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_equalizationCheckBox_clicked(bool checked)
{
    CameraOptions opt;
    opt.equalization = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_streamPortSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.portSendStream = ui->streamPortSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}



void CameraOptionsWindow::on_exposureHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.exposure = (double)value / 100.0;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->curExpLabel->setText(QString::number(opt.exposure));
}

void CameraOptionsWindow::on_pixelClockHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.pixelClock = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->curPixelClockLabel->setText(QString::number(opt.pixelClock));
}

void CameraOptionsWindow::on_frameRateHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.frameRate = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->curFrameRateLabel->setText(QString::number(opt.frameRate));
}

void CameraOptionsWindow::on_autoMinGainSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.autoExpParams.minGainCoeff = (double)ui->autoMinGainSpinBox->value(); //tmp
    cameraServer->sendParametersToCamera(opt, currentCamera);
}


void CameraOptionsWindow::on_autoMaxGainSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.autoExpParams.maxGainCoeff = (double)ui->autoMaxGainSpinBox->value();// tmp
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_autoPercentDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.autoExpParams.maxPercent = ui->autoPercentDoubleSpinBox->value() / 100.0;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_enableAutoExposureCheckBox_clicked(bool checked)
{
    CameraOptions opt;
    opt.internalAutoExposure = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->gainSlider->setDisabled(checked);
    ui->gainLineEdit->setDisabled(checked);
}

void CameraOptionsWindow::on_dayRadioButton_clicked(bool checked)
{
    if (checked)
    {
        CameraOptions opt;
        opt.autoExpParams.light = LightningParameter::Day;
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }

}

void CameraOptionsWindow::on_eveningRadioButton_clicked(bool checked)
{
    if (checked)
    {
        CameraOptions opt;
        opt.autoExpParams.light = LightningParameter::Evening;
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
}


void CameraOptionsWindow::on_saturationHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.saturation = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->saturationLineEdit->setText(QString::number(value));
}

void CameraOptionsWindow::on_saturationLineEdit_editingFinished()
{
    ui->saturationHorizontalSlider->setValue(ui->saturationLineEdit->text().toInt());
}

void CameraOptionsWindow::on_saturationRHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.rSaturation = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->saturationRLineEdit->setText(QString::number(value));
}

void CameraOptionsWindow::on_saturationRLineEdit_editingFinished()
{
    ui->saturationRHorizontalSlider->setValue(ui->saturationRLineEdit->text().toInt());
}

void CameraOptionsWindow::on_saturationGHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.gSaturation = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->saturationGLineEdit->setText(QString::number(value));
}

void CameraOptionsWindow::on_saturationGLineEdit_editingFinished()
{
    ui->saturationGHorizontalSlider->setValue(ui->saturationGLineEdit->text().toInt());
}

void CameraOptionsWindow::on_saturationBHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.bSaturation = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->saturationBLineEdit->setText(QString::number(value));
}

void CameraOptionsWindow::on_saturationBLineEdit_editingFinished()
{
    ui->saturationBHorizontalSlider->setValue(ui->saturationBLineEdit->text().toInt());
}


void CameraOptionsWindow::on_shadowCoefDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.shadowCoef = ui->shadowCoefDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}


void CameraOptionsWindow::on_shadowGaussWindowSizeSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.shadowGaussWindowSize = ui->shadowGaussWindowSizeSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_hueHorizontalSlider_valueChanged(int value)
{
    CameraOptions opt;
    opt.hue = value;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    ui->hueLineEdit->setText(QString::number(value));
}


void CameraOptionsWindow::on_hueLineEdit_editingFinished()
{
    ui->hueHorizontalSlider->setValue(ui->hueLineEdit->text().toInt());
}

void CameraOptionsWindow::on_autoLeftBorderDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.autoExpParams.minRelCoef = ui->autoLeftBorderDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_autoRightBorderDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.autoExpParams.maxRelCoef = ui->autoRightBorderDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}




void CameraOptionsWindow::on_shadowThresholdDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.shadowThreshold = ui->shadowThresholdDoubleSpinBox->value();
    qDebug() << "send shadow threshold";
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_chooseWBROIToolButton_clicked()
{
    QRect rect = roiScene->getROI();
    if (!rect.isEmpty())
    {
        ui->wbLabel->setStyleSheet("QLabel { background-color : green; }");
        CameraOptions opt;
        opt.wbRect = Rect(rect.x(), rect.y(), rect.width(), rect.height());
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
    roiScene->clearROI();
}

void CameraOptionsWindow::on_showWBROIToolButton_clicked()
{
    auto params = cameraServer->getLastCurrentParameters(currentCamera);
    if (!params.wbRect.empty())
    {
        roiScene->clearROI();
        QRect r = QRect(params.wbRect.x, params.wbRect.y, params.wbRect.width, params.wbRect.height);
        roiScene->drawROI(r);
    }


}

void CameraOptionsWindow::on_chooseMainROIPushButton_clicked()
{
    QRect rect = roiScene->getROI();
    if (!rect.isEmpty())
    {
        ui->mainROILabel->setStyleSheet("QLabel { background-color : green; }");
        CameraOptions opt;
        opt.recRois.mainSearchRect = Rect(rect.x(), rect.y(), rect.width(), rect.height());
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }

    roiScene->clearROI();
}

void CameraOptionsWindow::on_showMainROIPushButton_clicked()
{
    auto params = cameraServer->getLastCurrentParameters(currentCamera);
    if (!params.recRois.mainSearchRect.empty())
    {
        roiScene->clearROI();
        QRect r = QRect(params.recRois.mainSearchRect.x, params.recRois.mainSearchRect.y, params.recRois.mainSearchRect.width, params.recRois.mainSearchRect.height);
        roiScene->drawROI(r);
    }
}

void CameraOptionsWindow::on_chooseFirstROIPushButton_clicked()
{
    QRect rect = roiScene->getROI();
    if (!rect.isEmpty())
    {
        ui->firstROILabel->setStyleSheet("QLabel { background-color : green; }");
        CameraOptions opt;
        opt.recRois.trackFirstRect = Rect(rect.x(), rect.y(), rect.width(), rect.height());
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
    roiScene->clearROI();
}

void CameraOptionsWindow::on_showFirstROIPushButton_clicked()
{
    auto params = cameraServer->getLastCurrentParameters(currentCamera);
    if (!params.recRois.trackFirstRect.empty())
    {
        roiScene->clearROI();
        QRect r = QRect(params.recRois.trackFirstRect.x, params.recRois.trackFirstRect.y, params.recRois.trackFirstRect.width, params.recRois.trackFirstRect.height);
        roiScene->drawROI(r);
    }
}

void CameraOptionsWindow::on_chooseSecondROIPushButton_clicked()
{
    QRect rect = roiScene->getROI();
    if (!rect.isEmpty())
    {
        ui->secondROILabel->setStyleSheet("QLabel { background-color : green; }");
        CameraOptions opt;
        opt.recRois.trackSecondRect = Rect(rect.x(), rect.y(), rect.width(), rect.height());
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
}

void CameraOptionsWindow::on_showSecondROIPushButton_clicked()
{
    auto params = cameraServer->getLastCurrentParameters(currentCamera);
    if (!params.recRois.trackSecondRect.empty())
    {
        roiScene->clearROI();
        QRect r = QRect(params.recRois.trackSecondRect.x, params.recRois.trackSecondRect.y, params.recRois.trackSecondRect.width, params.recRois.trackSecondRect.height);
        roiScene->drawROI(r);
    }
}



void CameraOptionsWindow::on_rawTcpStreamCheckBox_clicked(bool checked)
{
    CameraOptions opt;
    opt.rawFrame = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_durationSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.videoDuration = ui->durationSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);

}

//tmp


void CameraOptionsWindow::on_triggerModeGroupBox_clicked(bool checked)
{
    CameraOptions opt;

    if (checked)
    {
        if (ui->highLowRadioButton->isChecked())
            opt.triggerMode = IS_SET_TRIGGER_HI_LO; // tmp
        else
            opt.triggerMode = IS_SET_TRIGGER_LO_HI;
    }
    else
    {
        opt.triggerMode = IS_SET_TRIGGER_OFF;
    }
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_highLowRadioButton_clicked(bool checked)
{
    if (checked)
    {
        CameraOptions opt;
        opt.triggerMode = IS_SET_TRIGGER_HI_LO;
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
}

void CameraOptionsWindow::on_lowHighRadioButton_clicked(bool checked)
{
    if (checked)
    {
        CameraOptions opt;
        opt.triggerMode = IS_SET_TRIGGER_LO_HI;
        cameraServer->sendParametersToCamera(opt, currentCamera);
    }
}

void CameraOptionsWindow::on_turnOnAutoVideoSaveRadioButton_clicked(bool checked)
{
    if (checked)
    {
        if (!ui->directoryLineEdit->text().isEmpty())
        {
            cameraServer->createVideoTimer(ui->timerSpinBox->value(), ui->durationSpinBox->value(), ui->fromTimeEdit->time(),
                                           ui->toTimeEdit->time(),ui->streamPortSpinBox->value(),
                                           currentCamera, ui->compressCheckBox->isChecked(), ui->directoryLineEdit->text());
        }
        else
        {
            QMessageBox::warning(this, "Внимание", "Не указана директория!");
        }
    }

}

void CameraOptionsWindow::on_turnOffAutoVideoSaveRadioButton_clicked(bool checked)
{
    if (checked)
    {
        cameraServer->stopVideoTimer(currentCamera);
    }

}

void CameraOptionsWindow::on_chooseDirectoryToolButton_clicked()
{
    static QString dir;
    dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                            "dir",
                                            QFileDialog::ShowDirsOnly
                                            | QFileDialog::DontResolveSymlinks);
    ui->directoryLineEdit->setText(dir);
}

void CameraOptionsWindow::on_corrDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.corrCoef = ui->corrDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}


void CameraOptionsWindow::on_searchAreaSizeSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.searchAreaSize = ui->searchAreaSizeSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_minSkoTemplateDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.minSkoOnTemplate = ui->minSkoTemplateDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_maxAngleBetwDirDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.maxAngleBetwDirections = ui->maxAngleBetwDirDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_minSpeedDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.minSpeed = ui->minSpeedDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_cannyMinSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.cannyThresMin = ui->cannyMinSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_cannyMaxSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.cannyThresMax = ui->cannyMaxSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_minAreaSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.minArea = ui->minAreaSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_maxAreaSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.maxArea = ui->maxAreaSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_circularityDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.circularityCoeff = ui->circularityDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_enableRecognizeCheckBox_toggled(bool checked)
{
    CameraOptions opt;
    opt.ballRecognizeFlag = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
    if (!checked)
    {
        cameraServer->clearRecognizeData();
    }
}

void CameraOptionsWindow::on_ballRecognizeStepSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.ballRecognizeStep = ui->ballRecognizeStepSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_debugRecFlagCheckBox_clicked(bool checked)
{
    CameraOptions opt;
    opt.debugRecFlag = checked;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_checkSynchroPushButton_clicked()
{
    if (!ui->directoryLineEdit->text().isEmpty())
    {
        cameraServer->syncVideo(ui->directoryLineEdit->text(), ui->compressCheckBox->isChecked());
    }
    else
    {
        QMessageBox::warning(this, "Внимание", "Не указан путь!");
    }

}

void CameraOptionsWindow::on_checkSyncPushButton_clicked()
{
    if (sync->isConnected())
    {
        sync->blockSignals(true);
        if (sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst, 1))
        {
            if (sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x1))
            {
                QThread::msleep(100);
                cameraServer->checkSync();
                QTimer::singleShot(500, this, [this]()
                {
                    sync->sendCommand(Synchronizer::SynchronizerProtocol::StopSync, 0x1);
                    sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst, 60);
                    sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x1);
                }
                );

            }
        }
        sync->blockSignals(false);
    }
    else
    {
        QMessageBox::warning(this, "Внимание", "Включите блок синхронизации!");
    }

}


void CameraOptionsWindow::on_restartCameraPushButton_clicked()
{
    cameraServer->restartCamera(currentCamera);
}

void CameraOptionsWindow::on_showVideoCheckBox_clicked(bool checked)
{
    cameraServer->setShowVideo(currentCamera, checked);
}

void CameraOptionsWindow::on_saveVideoCheckBox_clicked(bool checked)
{
    cameraServer->setWriteVideo(currentCamera, checked);
}

void CameraOptionsWindow::on_syncFramePushButton_clicked()
{
    cameraServer->syncFrame();
}

void CameraOptionsWindow::on_chooseReferenceFramePushButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/home",
                                                    tr(""));
    if (!fileName.isEmpty())
    {
        ui->referenceFrameLineEdit->setText(fileName);
    }
}



void CameraOptionsWindow::on_showPushButton_clicked()
{
    if (ui->referenceFrameLineEdit->text().isEmpty())
    {
        QMessageBox::warning(this, "Внимание", "Не указан один из путей");
    }
    else
    {
        Mat frame = cameraServer->getLastCameraFrame(currentCamera);
        if (frame.cols == 0 || frame.rows == 0)
        {
            QMessageBox::warning(this, "Внимание", "Запросите кадр с камеры!");
        }
        else
        {
            Mat fullFrame = Mat(Size(1936, 1216), CV_8UC3);
            frame.copyTo(fullFrame(Rect(0, 0, 1920, 1080)));
            Mat img = calibAdjustHelper.createCalibrateImage(ui->referenceFrameLineEdit->text(),
                                                             ui->streamPortSpinBox->value(),
                                                             frame, ui->patterWidthSpinBox->value(),
                                                             ui->searchWidthSpinBox->value(),
                                                             ui->searchHeightSpinBox->value());
            emit correlatedImageReady(img, currentCamera);
        }
    }
}

void CameraOptionsWindow::on_calibratePushButton_clicked()
{
    QVector <qint32> excludePoints;
    QStringList excludePointsStr = ui->excludePointsLineEdit->text().split(",");
    if (!excludePointsStr.isEmpty())
    {
        for (const auto& i : excludePointsStr)
        {
            excludePoints.append(i.toInt());
        }
    }
    Calibration::ExteriorOr eOr, nEOr;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cam, nCam;

    CalibrationAdjustHelper::readCurrentCalibrationParameters(ui->streamPortSpinBox->value(),
                                                              "calibrate/", eOr, cam, true);
    calibAdjustHelper.recalibrate(excludePoints, eOr,  cam, nEOr, nCam, true);
}

void CameraOptionsWindow::on_debounceEnableGroupBox_toggled(bool arg1)
{
    CameraOptions opt;
    opt.debounceEnable = arg1;
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_debounceSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.debounceValue = ui->debounceSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

void CameraOptionsWindow::on_saveCameraSettingsPushButton_clicked()
{
    cameraServer->saveCameraSettings(currentCamera);
}

void CameraOptionsWindow::on_skoCountDoubleSpinBox_editingFinished()
{
    CameraOptions opt;
    opt.recParams.skoCoef = ui->skoCountDoubleSpinBox->value();
    cameraServer->sendParametersToCamera(opt, currentCamera);
}

\

void CameraOptionsWindow::on_fullPictureCheckBox_clicked(bool checked)
{
    auto& av = ApproximationVisualizer::instance();
    av.createfullPicture(checked);
}

void CameraOptionsWindow::on_shortPictureCheckBox_clicked(bool checked)
{
    auto& av = ApproximationVisualizer::instance();
    av.createshortPicture(checked);
}

void CameraOptionsWindow::on_turnOnSyncFrameEveryCheckBox_toggled(bool checked)
{
    ApproximationVisualizer::instance().clearCalibGraphs();
    cameraServer->enableAutoCalibrate(checked, ui->syncFrameEverySpinBox->value(),
                                      ui->updateEVOAutoCalibrateCheckBox->isChecked(),
                                      flag);
}

void CameraOptionsWindow::loadSettings()
{
    QSettings settings;
    ui->turnOnSyncFrameEveryCheckBox->setChecked(settings.value("sync_frame_flag", false).toBool());
    ui->syncFrameEverySpinBox->setValue(settings.value("sync_frame_every", 10).toInt());

    flag = (CompareFlag)settings.value("compare_flag", CompareFlag::None).toInt();
    if (flag == None)
    {
        ui->compareNeibRadioButton->setChecked(true);
    }
    else if (flag == Reference)
    {
        ui->compareReferenceRadioButton->setChecked(true);
    }
    else
    {
        ui->compareCurrentRadioButton->setChecked(true);
    }
    cameraServer->enableAutoCalibrate(ui->turnOnSyncFrameEveryCheckBox->isChecked(),
                                      ui->syncFrameEverySpinBox->value(),
                                      ui->updateEVOAutoCalibrateCheckBox->isChecked(),
                                      flag);
    ApproximationVisualizer& av = ApproximationVisualizer::instance();
    av.createfullPicture(settings.value("big_picture", true).toBool());
    ui->fullPictureCheckBox->setChecked(av.isfullPicture());
    av.createshortPicture(settings.value("small_picture", true).toBool());
    ui->shortPictureCheckBox->setChecked(av.isshortPicture());
    av.setUpdateGraph(settings.value("update_graph", true).toBool());
    ui->updateGraphsCheckBox->setChecked(av.getUpdateGraph());
    ui->updateEVOAutoCalibrateCheckBox->setChecked(settings.value("update_evo_autocalibrate", false).toBool());

}

void CameraOptionsWindow::saveSettings()
{
    QSettings settings;
    settings.setValue("sync_frame_flag", ui->turnOnSyncFrameEveryCheckBox->isChecked());
    settings.setValue("sync_frame_every", ui->syncFrameEverySpinBox->value());
    ApproximationVisualizer& av = ApproximationVisualizer::instance();

    settings.setValue("big_picture", av.isfullPicture());
    settings.setValue("small_picture", av.isshortPicture());
    settings.setValue("update_graph",  av.getUpdateGraph());
    settings.setValue("compare_flag", flag);
    settings.setValue("update_evo_autocalibrate", ui->updateEVOAutoCalibrateCheckBox->isChecked());

}

void CameraOptionsWindow::on_updateGraphsCheckBox_toggled(bool checked)
{
    auto& av = ApproximationVisualizer::instance();
    av.setUpdateGraph(checked);
}

void CameraOptionsWindow::on_turnOnRecogAllCamerasCheckBox_toggled(bool checked)
{
    cameraServer->enableRecognitionAll(checked);
}



void CameraOptionsWindow::on_compareCurrentRadioButton_toggled(bool checked)
{
    flag = Current;
}

void CameraOptionsWindow::on_compareReferenceRadioButton_toggled(bool checked)
{
    flag = Reference;
}

void CameraOptionsWindow::on_compareNeibRadioButton_toggled(bool checked)
{
    flag = None;
}
