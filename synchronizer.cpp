#include "synchronizer.h"

Synchronizer::Synchronizer(qint32 port, QObject *parent) : QObject(parent)
{
    if (port != -1)
    {
        connect(port);
    }
}

bool Synchronizer::connect(qint32 port)
{
    QString error;
    if (port != -1 && QSerialPortInfo::availablePorts().size() > 0)
    {

        syncHard.setPort(QSerialPortInfo::availablePorts().at(port));
        syncHard.setBaudRate(9600);
        if (syncHard.open(QIODevice::ReadWrite))
        {
            qDebug() << "synchronizer connected";
            return true;
        }
        else
        {
            error = "Can't open connection on current port " + syncHard.errorString();
            emit errorOccured(error);
            qDebug() << error;
        }
    }
    else
    {
        error = "Invalid port or com ports not existss";
        qDebug() << error;
        emit errorOccured(error);
    }
    return false;
}

Synchronizer::~Synchronizer()
{
    disconnect();
}

bool Synchronizer::sendCommand(Synchronizer::SynchronizerProtocol command, quint8 value)
{

    QByteArray array;
    array.append(command);
    array.append(value);
    qDebug() << "write" << array.toHex();
    if (syncHard.write(array) == 2)
    {
        QByteArray rArray;
        while (syncHard.waitForReadyRead(500))
        {
            while (syncHard.bytesAvailable() > 0)
            {
                rArray.append(syncHard.readAll());
            }

        }
        qDebug() << "read" << rArray.toHex();
        if (rArray.isEmpty())
        {
            emit errorOccured("Не удалось получить ответное слово");
        }
        else
        {
            emit changed();
            return true;
        }
    }
    else
    {
        emit errorOccured(syncHard.errorString());
    }
    return false;
}

