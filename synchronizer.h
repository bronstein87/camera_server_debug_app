#ifndef SYNCHRONIZER_H
#define SYNCHRONIZER_H

#include <QObject>
#include <QtSerialPort/QSerialPort>
#include <QSerialPortInfo>
#include <QDebug>
class Synchronizer : public QObject
{
    Q_OBJECT
public:
    enum SynchronizerProtocol
    {
        CheckInterface = 0x00,
        TurnOn = 0x01,
        TurnOff = 0x02,
        StartSync = 0x03,
        StopSync = 0x04,
        SetFrameRateFirst = 0x0A,
        SetFrameRateSecond = 0x0B
    };

    explicit Synchronizer(qint32 port = -1, QObject *parent = 0);

    bool connect(qint32 port);

    bool isConnected() {return syncHard.isOpen();}

    void disconnect() {syncHard.close();}

    ~Synchronizer();

    bool sendCommand(SynchronizerProtocol command, quint8 value);
signals:
    void errorOccured(const QString& str);

    void changed();

private:
    QSerialPort syncHard;
};

#endif // SYNCHRONIZER_H
