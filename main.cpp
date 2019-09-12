#include "mainwindow.h"
#include <QApplication>
//#include <gst/gst.h>
#include <QSettings>

int main(int argc, char *argv[])
{
    //gst_init (&argc, &argv);
    QCoreApplication::setOrganizationName("IKI RAN");
    QCoreApplication::setApplicationName("BaseballCameraServer");
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}




