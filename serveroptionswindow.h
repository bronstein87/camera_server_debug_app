#ifndef SERVEROPTIONSWINDOW_H
#define SERVEROPTIONSWINDOW_H

#include <QWidget>
#include <synchronizer.h>
#include <QString>
#include <QSettings>


namespace Ui {
class ServerOptionsWindow;
}

class ServerOptionsWindow : public QWidget
{
    Q_OBJECT

public:
    explicit ServerOptionsWindow(Synchronizer* snc, QWidget *parent = 0);
    ~ServerOptionsWindow();

private slots:
    void on_turnOnSynchronizerPushButton_clicked();

    void on_synchroFirstCheckBox_clicked(bool checked);

    void on_synchroSecondCheckBox_clicked(bool checked);

    void on_turnOnFirstCheckBox_clicked(bool checked);

    void on_turnOnSecondCheckBox_clicked(bool checked);

    void on_disconnectPushButton_clicked();

    void on_allTurnOnCheckBox_clicked(bool checked);

    void on_frameRateFirstSpinBox_editingFinished();

    void on_frameRateSecondSpinBox_editingFinished();

private:

    void loadSettings();
    void saveSettings();

    Ui::ServerOptionsWindow *ui;
    Synchronizer* sync = nullptr;
};

#endif // SERVEROPTIONSWINDOW_H
