#ifndef CAMERASERVER_H
#define CAMERASERVER_H

#include <QObject>
#include <QTcpServer>
#include <QTcpSocket>
#include <QString>
#include <QVector>
#include <QScopedPointer>
#include <QHostAddress>
#include <QNetworkInterface>
#include <QMap>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
//#include <opencv2/cudaimgproc.hpp>
//#include <opencv2/cudaarithm.hpp>
//#include <opencv2/cudacodec.hpp>
#include <QtConcurrent/QtConcurrent>
#include <QLinkedList>
#include <QMutex>
#include <QMutexLocker>
#include <QTextStream>
#include <QFuture>
#include <calibration.h>
#include <ballapproximator.h>
#include <QXmlStreamReader>
#include <calibrationadjusthelper.h>
#include <QPixmap>
#include <QGraphicsPixmapItem>
#include <QPen>
#include <QGraphicsView>
#include <QFont>
#include <QFontDatabase>
#include <approximationvisualizer.h>
#include <limits.h>



enum Camera57Protocol
{
    SendCameraParams = 0x50,
    RequestCameraVideo = 0x56,
    RequestCameraFrame = 0x46,
    ReadyToGetStream = 0x52,
    IsServerPrepareToGetStream = 0x51,
    GetThrowCoordinates = 0x44,
    GetHitCoordinates = 0x45,
    GetTestDataFromClient = 0x54,
    StartStream = 0x53,
    StopStream = 0x55,
    RestartCamera = 0x59,
    AskCurrentCameraParams = 0x57,
    SendCurrentFrame = 0x60,
    SendSaveParameters = 0x61,
    SendRecognizeVideo = 0x62,
    EndOfMessageSign = 0x23
};



enum LightningParameter
{
    Invalid,
    Day,
    Evening
};

enum AppendFrameRule
{
    JustLast,
    JustQueue,
    All
};


#pragma pack(push, 1)
struct RecROIs
{
    cv::Rect mainSearchRect;
    cv::Rect trackFirstRect;
    cv::Rect trackSecondRect;
    cv::Rect mainHitSearchRect;
};
#pragma pack(pop)


#pragma pack(push, 1)
struct AutoExposureParameters
{
    qint32 light = Invalid;
    double gain = -1;
    double exposure = 1;
    double minGainCoeff = -1;
    double maxGainCoeff= -1;
    double maxPercent = -1;
    double minRelCoef = -1;
    double maxRelCoef = -1;
};

#pragma pack(pop)



#pragma pack(push, 1)
struct RecognizeParameters
{
    double corrCoef = -1;
    double skoCoef = -1; // mean + coef * sko for threshold for substracted frame
    qint32 searchAreaSize = -1; // search area for already tracked ball (prev_center - searchAreaSize, prev_center + searchAreaSize)
    double minSkoOnTemplate = -1;
    double maxAngleBetwDirections = -1; //degrees
    double minSpeed = -1; // pixel
    qint32 cannyThresMin = -1;
    qint32 cannyThresMax = -1;
    qint32 maxArea = -1;
    qint32 minArea = -1;
    double circularityCoeff = -1;
    double dirX = -1;
    double dirY = -1;
};
#pragma pack(pop)

static const constexpr int defaultHSV = -256;
#pragma pack(push, 1)
struct CurrentCameraParams
{
    double exposure;
    double minExposure;
    double maxExposure;
    qint32 pixelClock;
    qint32 minPixelClock;
    qint32 maxPixelClock;
    double frameRate;
    double minFrameRate;
    double maxFrameRate;
    qint32 sendFrameRate;
    qint32 gain;
    RecROIs recRois;
    cv::Rect wbRect; // tmp
    AutoExposureParameters autoExpParams;
    qint32 rawFrame = -1;
    qint32 videoDuration;
    qint32 width;
    qint32 height;
    qint32 triggerMode = 0;
    qint32 pictureParamFlag = 0;
    qint32 whiteBalance = 0;
    qint32 equalization = 0;
    qint32 autoExposureIntenal = 0;
    qint32 focusing = -1;
    double gamma = -1;
    double sharp = -1;
    qint32 portSendStream = -1;
    qint32 IRCorrection = 0;
    qint32 saturation = defaultHSV;
    qint32 rSaturation = defaultHSV;
    qint32 gSaturation = defaultHSV;
    qint32 bSaturation = defaultHSV;
    qint32 hue = defaultHSV;
    double shadowCoef = -1;
    double shadowThreshold = -1;
    qint32 shadowGaussWindowSize = -1;
    RecognizeParameters recParams;
    qint32 ballRecognizeFlag = -1;
    qint32 ballRecognizeStep = -1;
    qint32 debugRecFlag = 0;
    qint32 debounceEnable = -1;
    qint32 debounceValue = -1;
    qint32 rotate = -1;
};
#pragma pack(pop)




constexpr const static qint32 maxNumberOfMeasures = 100;

constexpr const static qint32 measureDim = 4;

struct RecognizeData
{
    double data[maxNumberOfMeasures][measureDim];
    qint32 size;
    QTime startTime;
};

struct RecognizeVideoData
{
    QTime startTime;
    qint32 frameCount;
    QVector <quint64> times;
};

struct CameraStatus
{
    bool streamIsActive = false;
    QLinkedList <cv::Mat> frameBuffer;
    QByteArray buffer;
    CurrentCameraParams curParams;
    QTimer* timer = nullptr;
    quint64 syncTime;
    QTime machineSyncTime;
    bool showVideo = false;
    bool writeVideo = false;
    QVector <quint64> times;
    cv::Mat lastFrame;
    QVector <RecognizeData> recData;
    QVector <RecognizeVideoData> recVideoData;
    QVector <RecognizeData> hitData;
    double firstTimeBall;

};



#pragma pack(push, 1)
struct CameraOptions
{
    qint32 pictureParamFlag = -1;
    qint32 whiteBalance = -1;
    qint32 equalization = -1;
    qint32 internalAutoExposure = -1;
    qint32 focusing = -1;
    double exposure = -1.0;
    double frameRate = -1.0;
    double frameRateSend = -1;
    qint32 pixelClock = -1;
    double gamma = -1;
    qint32 gain = -1;
    double sharp = -1;
    qint32 AOIHeight = -1;
    qint32 AOIWidth = -1;
    qint32 portSendStream = -1;
    qint32 IRCorrection = -1;
    qint32 saturation = defaultHSV;
    qint32 rSaturation = defaultHSV;
    qint32 gSaturation = defaultHSV;
    qint32 bSaturation = defaultHSV;
    qint32 hue = defaultHSV;
    double shadowCoef = -1;
    double shadowThreshold = -1;
    qint32 shadowGaussWindowSize = -1;
    RecROIs recRois;
    cv::Rect wbRect; // tmp
    AutoExposureParameters autoExpParams;
    qint32 rawFrame = -1;
    qint32 videoDuration = -1;
    qint32 triggerMode = -1;
    RecognizeParameters recParams;
    qint32 ballRecognizeFlag = -1;
    qint32 ballRecognizeStep = -1;
    qint32 debugRecFlag = -1;
    qint32 debounceEnable = -1;
    qint32 debounceValue = -1;
    qint32 rotate = -1;
};
#pragma pack(pop)


class CameraServer : public QObject
{
    Q_OBJECT
public:
    explicit CameraServer(QObject *parent = nullptr);

    ~CameraServer();

    void startServer();

    void requestVideoFromCamera(qint32 port, qint32 duration, bool compress, QTcpSocket* socket, const QString& path);

    void requestFrameFromCamera(QTcpSocket* socket, qint32 port, const QString &path);

    void sendParametersToCamera(const CameraOptions& parameters, QTcpSocket* socket);

    void receiveRecognizeVideo(QTcpSocket* camera);

    void requestCurrentCameraParameters(QTcpSocket* socket);

    void updateParamsForAllCameras();

    void startStream(qint32 port, QTcpSocket* socket, bool start);

    void restartCamera(QTcpSocket* socket);

    void saveCameraSettings(QTcpSocket* camera);

    void checkRecognizeResults();

    void checkHitResults();

    cv::Mat getFrameFromStream(QTcpSocket* socket);

    cv::Mat getLastCameraFrame(QTcpSocket* camera);

     cv::Mat getLastCameraFrame(qint32 num);

    CurrentCameraParams& getLastCurrentParameters(QTcpSocket* socket) {return cameras[socket].curParams;}

    void createVideoTimer(qint32 interval, qint32 duration, QTime from, QTime to, qint32 port, QTcpSocket* socket, bool compress, const QString &path);

    void stopVideoTimer(QTcpSocket* socket);

    void syncVideo(const QString& dir, bool compress);

    void checkSync();

    void syncFrame(AppendFrameRule rule = All, bool dontSetStream = false);

    void setSyncAnyWay(bool f) {syncAnyWay = f;}

    void testApproximation(const QString& filePath, BallApproximator &approx, const QString &prefix = QString(), double coef = 0.04);

    void setWriteVideo(QTcpSocket* socket, bool f) {cameras[socket].writeVideo = f;}

    void setShowVideo(QTcpSocket* socket, bool f) {cameras[socket].showVideo = f;}

    bool getWriteVideo(QTcpSocket* socket) {return cameras[socket].writeVideo;}

    bool getShowVideo(QTcpSocket* socket) {return cameras[socket].showVideo;}

    void enableAutoCalibrate(bool flag, qint32 syncFrameEvery, bool save, CompareFlag cFlag);

    bool getAutoCalibrate() {return autoCalibrate;}

    void enableRecognitionAll(bool flag);

    void runRepeatStream();

    void clearRecognizeData()
    {
        for (auto& i : cameras)
        {
            i.recData.clear();
            i.recVideoData.clear();
        }
    }


signals:

    void readyMessageFromServer(const QString& msg, QTcpSocket* client = nullptr);

    void connectedCameraSocket(QTcpSocket* socket);

    void disconnectedCameraSocket(QTcpSocket* socket);

    void streamFromIsComing(QTcpSocket* socket);

    void frameFromStreamReady(QTcpSocket* socket);

    void finishedGetData(QTcpSocket* socket);

    void currentCameraParamsReady(QTcpSocket* camera, const CurrentCameraParams& params);

    void syncTimeGot();

    void resultPictureReady(QImage& image);

    void autoCalibrateImageReady(cv::Mat mat, QTcpSocket* socket);

private:

    struct WriteRawVideoData
    {
        WriteRawVideoData(QDataStream& _in, QTextStream& _ints, cv::VideoWriter& _video, qint32 _width, qint32 _height)
            :in(_in), ints(_ints), video(_video), width(_width), height(_height) {}
        QDataStream &in;
        QTextStream &ints;
        cv::VideoWriter &video;
        qint32 width;
        qint32 height;
    };

    bool handleApproximation(BallApproximator& approx, qint32 fNum, double points1[maxNumberOfMeasures][measureDim], qint32 size1,  qint32 sNum, double points2[maxNumberOfMeasures][measureDim], qint32 size2);

    void correctSyncTime();

    void startServerInternal(QTcpServer* server, qint32 port);

    void getFrameInternal(QTcpSocket* camera, qint32 port, const QString& path, bool sync, AppendFrameRule rule = All, bool dontSetStream = false);

    void getVideoInternal(QTcpSocket *camera, qint32 frameCount, qint32 port, bool compress, const QString& path);

    void handleNewConnection();

    void handleMessageFromCamera(QTcpSocket* camera);

    void appendFrameToBuffer(cv::Mat frame, QTcpSocket* camera, AppendFrameRule rule = All);

    void writeRawVideoStream(QTcpSocket *camera, QByteArray &buffer, cv::Mat matRaw, cv::Mat mat, qint32 &i, WriteRawVideoData& d);

    void fillEndOfMessage(QByteArray& array);

    void receiveRtspVideo(qint32 port, qint32 frameCount, const QString& path, QTcpSocket* camera, const QVector<quint64> &times = QVector<quint64>());

    qint32 correctCameraTime(quint64 &intTime, QTcpSocket* camera);

    QScopedPointer <QTcpServer> commandServer;
    QMap <QTcpSocket*, CameraStatus> cameras;
    QVector <char> avaliableCommands;
    QMutex streamMutex; // сделать mutex для каждой камеры?
    QMutex autoCalibrateMutex;
    QMutex repeatMutex;
    bool calibrate = false;
    QMap <qint32, QAtomicInteger <qint64>> syncTime;
    bool syncAnyWay = false;
   // QAtomicInteger <qint64> syncTime = 0;
    //QAtomicInteger <qint8> repeatCameraChosen = 0;
    bool autoCalibrate = false;
    bool updateAutoCalibrate = false;
    CompareFlag compareFlag;
    qint32 syncFrameEvery = 10;
    QTimer frameEveryTimer;
    qint32 waitFor;


    static constexpr const qint32 exposureToCameraIntTime = 10000;
    static constexpr const qint32 correctSyncTimeEvery = 100;
    static constexpr const qint32 machineTimeEpsilon = 100;
    static constexpr const qint32 syncTimeEpsilon = 5000;



};

#endif // CAMERASERVER_H
