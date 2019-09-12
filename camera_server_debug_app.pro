#-------------------------------------------------
#
# Project created by QtCreator 2019-05-30T09:18:02
#
#-------------------------------------------------

QT += core gui network serialport
QT += printsupport
QT += gui widgets
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = camera_server_debug_app
TEMPLATE = app
CONFIG += c++14

SOURCES += main.cpp\
        mainwindow.cpp \
    cameraserver.cpp \
    editframeqgraphicsview.cpp \
    zoomgraphicsview.cpp \
    cameraoptionswindow.cpp \
    synchronizer.cpp \
    serveroptionswindow.cpp \
    calibration.cpp \
    ballapproximator.cpp \
    calibrationadjusthelper.cpp \
    mathfunc.cpp \
    approximationvisualizer.cpp \
    qcustomplot/qcustomplot.cpp \
    qcustomplot/cplotter.cpp \
    qcustomplot/cxyplotter.cpp

HEADERS  += mainwindow.h \
    cameraserver.h \
    editframeqgraphicsview.h \
    zoomgraphicsview.h \
    cameraoptionswindow.h \
    synchronizer.h \
    serveroptionswindow.h \
    calibration.h \
    ballapproximator.h \
    calibrationadjusthelper.h \
    mathfunc.h \
    approximationvisualizer.h \
    qcustomplot/qcustomplot.h \
    qcustomplot/cplotter.h \
    qcustomplot/cxyplotter.h \
    goodcolors.h

FORMS    += mainwindow.ui \
    cameraoptionswindow.ui \
    serveroptionswindow.ui


#INCLUDEPATH += $$PWD/../../../usr/local/include
#DEPENDPATH += $$PWD/../../../usr/local/include



#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_core

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaarithm

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaimgproc

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_shape

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_highgui

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaimgproc

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_imgcodecs

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_features2d

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudacodec

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudalegacy

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_videoio

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudafilters

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_video

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_imgproc

#unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_features2d



#INCLUDEPATH += /usr/lib/aarch64-linux-gnu/gstreamer-1.0/include
#INCLUDEPATH += /usr/include/glib-2.0
#INCLUDEPATH += /usr/lib/aarch64-linux-gnu/glib-2.0/include
#INCLUDEPATH += /usr/include/gstreamer-1.0sa

#unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/aarch64-linux-gnu/ -lgstreamer-1.0

#unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/aarch64-linux-gnu/ -lgobject-2.0

#unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/aarch64-linux-gnu/ -lglib-2.0

#unix:!macx: LIBS += -L$$PWD/../../../../usr/lib/aarch64-linux-gnu/ -lgio-2.0




#win32: LIBS += -LD:/gstreamer/1.0/x86_64/lib/ -lgstreamer-1.0

#INCLUDEPATH += D:/gstreamer/1.0/x86_64/include
#DEPENDPATH += D:/gstreamer/1.0/x86_64/include
#INCLUDEPATH += D:/gstreamer/1.0/x86_64/include/gstreamer-1.0
#DEPENDPATH += D:/gstreamer/1.0/x86_64/include/gstreamer-1.0

#win32: LIBS += -LD:/gstreamer/1.0/x86_64/lib/ -lglib-2.0

#INCLUDEPATH += D:/gstreamer/1.0/x86_64/include/glib-2.0
#DEPENDPATH += D:/gstreamer/1.0/x86_64/include/glib-2.0

#INCLUDEPATH += $$PWD/../../../gstreamer/1.0/x86_64/lib/glib-2.0/include
#DEPENDPATH += $$PWD/../../../gstreamer/1.0/x86_64/lib/glib-2.0/include


#win32: LIBS += -LD:/gstreamer/1.0/x86_64/lib/ -lgobject-2.0

#INCLUDEPATH += D:/gstreamer/1.0/x86_64/lib/glib-2.0/include
#DEPENDPATH += D:/gstreamer/1.0/x86_64/lib/glib-2.0/include







RESOURCES += \
    resources.qrc














win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_core346
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_core346d


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_imgproc346
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_imgproc346d



INCLUDEPATH += $$PWD/../../build_opencv/install/include
DEPENDPATH += $$PWD/../../build_opencv/install/include



win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_video346
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_video346d


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_videoio346
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_videoio346d


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_imgcodecs346
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../build_opencv/install/x64/vc15/lib/ -lopencv_imgcodecs346d


