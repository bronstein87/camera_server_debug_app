#-------------------------------------------------
#
# Project created by QtCreator 2019-05-30T09:18:02
#
#-------------------------------------------------

QT       += core gui network serialport

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
    mathfunc.cpp

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
    mathfunc.h

FORMS    += mainwindow.ui \
    cameraoptionswindow.ui \
    serveroptionswindow.ui


INCLUDEPATH += $$PWD/../../../usr/local/include
DEPENDPATH += $$PWD/../../../usr/local/include

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_core

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaarithm

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaimgproc

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_shape

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_highgui

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudaimgproc

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_imgcodecs

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_features2d

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudacodec

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudalegacy

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_videoio

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_cudafilters

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_video

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_imgproc

unix:!macx: LIBS += -L$$PWD/../../../usr/local/lib/ -lopencv_features2d
