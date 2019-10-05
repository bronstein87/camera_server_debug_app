#include "serveroptionswindow.h"
#include "ui_serveroptionswindow.h"

ServerOptionsWindow::ServerOptionsWindow(Synchronizer* snc, QWidget *parent) :
    sync(snc), QWidget(parent),
    ui(new Ui::ServerOptionsWindow)
{
    ui->setupUi(this);

    setWindowFlags(Qt::WindowStaysOnTopHint);
    ui->statusLabel->setText("Выключен");
    ui->statusLabel->setStyleSheet("QLabel { background-color : red; }");
    connect(sync, &Synchronizer::errorOccured, [this](auto str)
    {
        ui->mainMonitorTextEdit->append(str);
    });
    ui->mainMonitorTextEdit->append("Доступные порты:");
    for (const auto& i : QSerialPortInfo::availablePorts())
    {
        QString str ("%1 %2 %3 %4");
        str = str.arg(i.description())
                .arg(i.manufacturer())
                .arg(i.portName())
                .arg(i.productIdentifier());
        ui->mainMonitorTextEdit->append(str);
    }
    loadSettings();
}

ServerOptionsWindow::~ServerOptionsWindow()
{
    saveSettings();
    delete ui;
}

void ServerOptionsWindow::on_turnOnSynchronizerPushButton_clicked()
{
    if (sync->connect(ui->portSpinBox->value()))
    {
        ui->statusLabel->setText("Подключено");
        ui->statusLabel->setStyleSheet("QLabel { background-color : green; }");
    }

}

void ServerOptionsWindow::on_synchroFirstCheckBox_clicked(bool checked)
{
    static bool previous = false;
    if (checked && !previous)
    {

        if (sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst,
                              ui->frameRateFirstSpinBox->value()))
        {
            if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x1))
            {
                ui->synchroFirstCheckBox->setChecked(false);
            }
            else
            {
                qDebug() << "synchro first" << ui->frameRateFirstSpinBox->value();
                previous = checked;
            }
        }
        else
        {

            ui->synchroFirstCheckBox->setChecked(false);
        }

    }
    else if (!checked && previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::StopSync, 0x1))
        {
            ui->synchroFirstCheckBox->setChecked(true);
        }
        else
        {
            previous = checked;
        }
    }
}

void ServerOptionsWindow::on_synchroSecondCheckBox_clicked(bool checked)
{
    static bool previous = false;
    if (checked && !previous)
    {
        if (sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateSecond,
                              ui->frameRateSecondSpinBox->value()))
        {
            if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::StartSync, 0x2))
            {
                ui->synchroSecondCheckBox->setChecked(false);
            }
            else
            {
                qDebug() << "synchro second" << ui->frameRateSecondSpinBox->value();
                previous = checked;
            }
        }
        else
        {
            ui->synchroFirstCheckBox->setChecked(false);
        }

    }
    else if (!checked && previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::StopSync, 0x2))
        {
            ui->synchroSecondCheckBox->setChecked(true);
        }
        else
        {
            previous = checked;
        }
    }
}

void ServerOptionsWindow::on_turnOnFirstCheckBox_clicked(bool checked)
{
    static bool previous = false;
    if (checked && !previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOn, 0x1))
        {
            ui->turnOnFirstCheckBox->setChecked(false);
        }
        else
        {
            qDebug() << "turn on first";
            previous = checked;
        }
    }
    else if (!checked && previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOff, 0x1))
        {
            ui->turnOnFirstCheckBox->setChecked(true);
        }
        else
        {
            previous = checked;
        }
    }
}

void ServerOptionsWindow::on_turnOnSecondCheckBox_clicked(bool checked)
{
    static bool previous = false;
    if (checked && !previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOn, 0x2))
        {
            ui->turnOnSecondCheckBox->setChecked(false);
        }
        else
        {
            qDebug() << "turn on second";
            previous = checked;
        }
    }
    else if (!checked && previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOff, 0x2))
        {
            ui->turnOnSecondCheckBox->setChecked(true);
        }
        else
        {
            previous = checked;
        }
    }
}

void ServerOptionsWindow::on_disconnectPushButton_clicked()
{
    sync->disconnect();
    ui->statusLabel->setStyleSheet("QLabel { background-color : red; }");
}



void ServerOptionsWindow::on_allTurnOnCheckBox_clicked(bool checked)
{
    static bool previous = false;
    if (checked && !previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOn, 0x3))
        {
            ui->turnOnSecondCheckBox->setChecked(false);
        }
        else
        {
            previous = checked;
        }
    }
    else if (!checked && previous)
    {
        if (!sync->sendCommand(Synchronizer::SynchronizerProtocol::TurnOff, 0x3))
        {
            ui->turnOnSecondCheckBox->setChecked(true);
        }
        else
        {
            previous = checked;
        }
    }
}

void ServerOptionsWindow::on_frameRateFirstSpinBox_editingFinished()
{
    sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateFirst,
                      ui->frameRateFirstSpinBox->value());
}

void ServerOptionsWindow::on_frameRateSecondSpinBox_editingFinished()
{
    sync->sendCommand(Synchronizer::SynchronizerProtocol::SetFrameRateSecond,
                      ui->frameRateSecondSpinBox->value());
}

void ServerOptionsWindow::loadSettings()
{
    QSettings settings;

    ui->turnOnActivateCheckBox->setChecked(settings.value("enable_sync", false).toBool());
    if (ui->turnOnActivateCheckBox->isChecked())
    {
        ui->portSpinBox->setValue(settings.value("com_port_num", 0).toInt());
        on_turnOnSynchronizerPushButton_clicked();
    }
}

void ServerOptionsWindow::saveSettings()
{
    QSettings settings;
    settings.setValue("enable_sync", ui->turnOnActivateCheckBox->isChecked());
    settings.setValue("com_port_num", ui->portSpinBox->value());
}
