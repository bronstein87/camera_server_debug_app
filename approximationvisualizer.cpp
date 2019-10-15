#include "approximationvisualizer.h"
using namespace cv;
ApproximationVisualizer::ApproximationVisualizer(QObject *parent) : QObject(parent),
    scenePicture (new QGraphicsScene()),
    fullPixmap("://resources/picture.png"), shortPixmap("://resources/picture3.png")
{
    //  batterSearchZones.insert(3850, cv::Rect(930, 866, 1250 - 930, 1050 - 866)); //tmp
    // batterSearchZones.insert(4510, cv::Rect(640, 888, 1035 - 640, 1045 - 888));
    batterSearchZones.insert(3850, cv::Rect(908, 644, 1244 - 908, 1035 - 644)); //tmp
    batterSearchZones.insert(4510, cv::Rect(671, 647, 966 - 671, 1044 - 647));
    scenarioData.scaleRects.insert(3850, cv::Rect(570, 230, 1745 - 570, 1032 - 230));
    scenarioData.scaleRects.insert(4510, cv::Rect(179, 188, 1383 - 179, 1040 - 188));
    auto img = shortPixmap.toImage();
    Mat m(img.height(), img.width(), CV_8UC4, (void*)img.constBits());
    cvtColor(m, m, CV_RGBA2RGB);
    Mat scaled;
    cv::resize(m, scaled, Size(scenarioWidth, scenarioHeight));
    scenarioData.initPicture = scaled;
}

QVector<QImage> ApproximationVisualizer::createApproxVisualisation(BallApproximator& approx, QString info)
{
    QVector <QImage> images;
    if (fullPicture)
    {
        images.append(makeFullPicture(approx));
    }
    if (shortPicture)
    {
        images.append(makeShortPicture(approx, info));
    }
    return images;
}

void ApproximationVisualizer::fillSkippedTime(QVector<double>& errorsFirst, QVector<double>& timeFirstInit, QVector<double>& timeFirst)
{
    qint32 from = 0;
    for (qint32 i = 0; i < timeFirstInit.size(); ++i)
    {
        timeFirstInit[i] = std::abs(timeFirstInit[i]);
        for (qint32 j = from; j < timeFirst.size(); ++j)
        {
            if (qFuzzyCompare(std::abs(timeFirstInit[i]), timeFirst[j]))
            {
                from = j + 1 < timeFirst.size() ? j + 1 : j;
                break;
            }
            if (timeFirstInit[i] < std::abs(timeFirst[j]))
            {
                from = j + 1;
                timeFirst.insert(j, timeFirstInit[i]);
                errorsFirst.insert(j, errorsFirst[j]);
                break;
            }
            if (timeFirstInit[i] > std::abs(timeFirst[j]))
            {
                from = j + 1;
                timeFirst.insert(j + 1, timeFirstInit[i]);
                errorsFirst.insert(j + 1, errorsFirst[j]);
                break;
            }
        }
    }
}

void ApproximationVisualizer::clear(bool all)
{
    if (all)
    {
        for (auto & i : plots)
        {
            i->clearGraphs();

        }
    }
    else
    {
        for (qint32 i = 0; i < TimeGraph + 1; ++i)
        {
            plots[i]->clearGraphs();
        }
    }


}

void ApproximationVisualizer::plotCamera(BallApproximator& approx, qint32 firstNumber, qint32 secondNumber, PlotCameraData& plotData, QDateTime dt)
{
    QVector <double> timeFirst = plotData.timeFirst;
    QVector <double> timeSecond = plotData.timeSecond;
    QVector <double> errorsFirst = approx.getFirstErrors();
    QVector <double> errorsSecond = approx.getSecondErrors();
    QVector <double> timeFirstDelta;
    QVector <double> timeSecondDelta;


    for (qint32 i = 1; i < errorsFirst.size(); ++i)
    {
        timeFirstDelta.append((timeFirst[i] - timeFirst [i - 1]) * 1000000);
    }
    timeFirstDelta.prepend(timeFirstDelta.first());
    for (qint32 i = 1; i < errorsSecond.size(); ++i)
    {
        timeSecondDelta.append((timeSecond[i] - timeSecond[i - 1]) * 1000000);
    }
    timeSecondDelta.prepend(timeSecondDelta.first());

    Measure m;
    m.error = errorsFirst;
    m.time = dt;
    m.timeDelta = timeFirstDelta;
    m.objCount = plotData.firstMesCount;
    measuresMap[firstNumber].append(m);

    m.error = errorsSecond;
    m.timeDelta = timeSecondDelta;
    m.objCount = plotData.secondMesCount;
    measuresMap[secondNumber].append(m);

    qint32 index = measuresMap[firstNumber].size() - 1;
    if (updateGraph)
    {
        clear(false);
        QVector <double> pointsFirst;
        for (qint32 i = 0; i < measuresMap[firstNumber][index].error.size(); ++i)
        {
            pointsFirst.append(i);
        }

        QVector <double> pointsSecond;
        for (qint32 i = 0; i < measuresMap[secondNumber][index].error.size(); ++i)
        {
            pointsSecond.append(i);
        }
        plotter.setPlot(plots[TimeGraph]);
        plotter.setTitle("Times");
        plotter.setGraphName(QString::number(firstNumber));
        plotter.setCustomColor(colorMap[firstNumber]);
        plotter.addDefaultGraph(measuresMap[firstNumber][index].timeDelta, pointsFirst);
        plotter.setGraphName(QString::number(secondNumber));
        plotter.setCustomColor(colorMap[secondNumber]);
        plotter.addDefaultGraph(measuresMap[secondNumber][index].timeDelta, pointsSecond);
        plots[TimeGraph]->saveJpg(QString("games_%1/graphs/times_%2.jpg")
                                  .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                  .arg(QTime::currentTime().toString("hh_mm_ss")), 400, 400);


        plotter.setPlot(plots[ErrorGraph]);
        plotter.setTitle("Errors");
        plotter.setGraphName(QString::number(firstNumber));
        plotter.setCustomColor(colorMap[firstNumber]);
        plotter.addDefaultGraph(measuresMap[firstNumber][index].error, pointsFirst);
        plotter.setGraphName(QString::number(secondNumber));
        plotter.setCustomColor(colorMap[secondNumber]);
        plotter.addDefaultGraph(measuresMap[secondNumber][index].error, pointsSecond);
        plots[ErrorGraph]->saveJpg(QString("games_%1/graphs/errors_%2.jpg")
                                   .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                   .arg(QTime::currentTime().toString("hh_mm_ss")), 400, 400);

        qint32 indexPlot = 0;
        for (auto& i : measuresMap.keys())
        {
            plotter.setPlot(plots[ObjCountGraph]);
            plotter.setTitle("Object count");
            plotter.setGraphName(QString::number(i));
            plotter.setCustomColor(colorMap[i]);
            plotter.addPoint(QDateTime::currentDateTime().toMSecsSinceEpoch(), measuresMap[i][index].objCount, PlotType::DATETIME, indexPlot);
            ++indexPlot;
        }
        plots[ObjCountGraph]->rescaleAxes();
        plots[ObjCountGraph]->replot();
        plots[ObjCountGraph]->saveJpg(QString("games_%1/graphs/obj_count_%2.jpg")
                                      .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                      .arg(QTime::currentTime().toString("hh_mm_ss")), 400, 400);
    }
    measureCountChanged(index);

}

void ApproximationVisualizer::plotCalibrate(QMap <qint32, Calibration::ExteriorOr>& map, CompareFlag compareFlag)
{

    QDateTime dt = QDateTime::currentDateTime();
    CalibData d;
    for (auto& i : map.keys())
    {
        d.dt = dt;
        d.eOr = map[i];
        calibMap[i].append(d);
    }
    QVector <QCustomPlot*> calibPlots {plots[CalibFirstGraph], plots[CalibSecondGraph]};
    qint32 indexPlot = 0;

    for (auto &i : calibMap.keys())
    {
        auto& v = calibMap[i];
        Calibration::ExteriorOr eo;
        Calibration::SpacecraftPlatform::CAMERA::CameraParams cam;
        if (compareFlag == Current)
        {
            CalibrationAdjustHelper::readCurrentCalibrationParameters(i, "calibrate/", eo, cam, true);
        }
        else if (compareFlag == Reference)
        {
            CalibrationAdjustHelper::readCurrentCalibrationParameters(i, "calibrate/", eo, cam, false);
        }
        if (!v.isEmpty())
        {
            qint32 index = v.size() - 1;
            plotter.setTitle(QString("Калибровка %1").arg(i));
            plotter.setPlot(calibPlots[indexPlot]);
            plotter.setGraphName("Omega");
            plotter.setCustomColor(QRgb(0x209fdf));
            double value;
            if (compareFlag != None)
            {
                value  = (v[index].eOr.OPK.Omega - eo.OPK.Omega) * BOKZMath::radToDegrees;
            }
            else
            {
                value  = index == 0 ? 0 : (v[index].eOr.OPK.Omega - v[index - 1].eOr.OPK.Omega) * BOKZMath::radToDegrees;
            }

            plotter.addPoint(v.last().dt.toMSecsSinceEpoch(), value, PlotType::DATETIME, 0);

            plotter.setPlot(calibPlots[indexPlot]);
            plotter.setGraphName("Phi");
            plotter.setCustomColor(QRgb(0x99ca53));
            if (compareFlag != None)
            {
                value = (v[index].eOr.OPK.Phi - eo.OPK.Phi) * BOKZMath::radToDegrees;
            }
            else
            {
                value = index == 0 ? 0 : (v[index].eOr.OPK.Phi - v[index - 1].eOr.OPK.Phi) * BOKZMath::radToDegrees;
            }
            plotter.addPoint(v.last().dt.toMSecsSinceEpoch(), value, PlotType::DATETIME, 1);

            //            plotter.setPlot(calibPlots[indexPlot]);
            //            plotter.setGraphName("Kappa");
            //            plotter.setCustomColor(QRgb(0xf6a625));
            //            if (compareWithStandart)
            //            {
            //                value = index == 0 ? 0 : (v[index].eOr.OPK.Kappa - eo.OPK.Kappa) * BOKZMath::radToDegrees;
            //            }
            //            else
            //            {
            //                value = index == 0 ? 0 : (v[index].eOr.OPK.Kappa - v[index - 1].eOr.OPK.Kappa) * BOKZMath::radToDegrees;
            //            }
            //            plotter.addPoint(v.last().dt.toMSecsSinceEpoch(), value, PlotType::DATETIME, 2);

            plotter.setPlot(calibPlots[indexPlot]);
            plotter.setGraphName("XYZ");
            plotter.setCustomColor(QRgb(0x6d5fd5));
            if (compareFlag != None)
            {

                double lenghtPrev =  sqrt(pow(eo.Point.X, 2)  + pow(eo.Point.Y, 2)
                                          + pow(eo.Point.Z, 2));
                double lenght =  sqrt(pow(v[index].eOr.Point.X, 2)  + pow(v[index].eOr.Point.Y, 2)
                                      + pow(v[index].eOr.Point.Z, 2));
                value = lenght - lenghtPrev;
            }
            else
            {
                if (index == 0)
                {
                    value = 0;
                }
                else
                {
                    double lenghtPrev =  sqrt(pow(v[index - 1].eOr.Point.X, 2)  + pow(v[index - 1].eOr.Point.Y, 2)
                            + pow(v[index - 1].eOr.Point.Z, 2));
                    double lenght =  sqrt(pow(v[index].eOr.Point.X, 2)  + pow(v[index].eOr.Point.Y, 2)
                                          + pow(v[index].eOr.Point.Z, 2));
                    value = lenght - lenghtPrev;
                }
            }

            plotter.addPoint(v.last().dt.toMSecsSinceEpoch(), value, PlotType::DATETIME, 2);
            calibPlots[indexPlot]->rescaleAxes();
            calibPlots[indexPlot]->replot();

            calibPlots[indexPlot]->saveJpg(QString("games_%1/graphs/obj_count_%2_%3.jpg")
                                           .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                           .arg(i)
                                           .arg(QTime::currentTime().toString("hh_mm_ss")), 400, 400);

            ++indexPlot;
        }
    }
}

void ApproximationVisualizer::runRepeatStream()
{
    if (!repeatThreadRun)
    {
        QtConcurrent::run(this, &ApproximationVisualizer::runRepeatStreamInternal);
    }
}

void ApproximationVisualizer::stopRepeatStream()
{
    repeatThreadRun = 0;
}

void ApproximationVisualizer::setPlots(const QVector<QCustomPlot *> &_plots)
{
    plots = _plots;
    for (qint32 i = ObjCountGraph; i <  CalibSecondGraph + 1; ++i)
    {
        plotSyncLines.append(new QCPItemStraightLine(plots[i]));
        plotSyncLines.last()->point1->setCoords(0, 0);
        plotSyncLines.last()->point2->setCoords(0, 0);
    }
}

void ApproximationVisualizer::initCamerasPlots(const QMap<qint32, QColor> &map)
{
    colorMap = map;
    for (auto& i : colorMap.keys())
    {
        measuresMap.insert(i, QVector <Measure> ());
        calibMap.insert(i, QVector <CalibData>());
    }
}


void ApproximationVisualizer::showMeasure(qint32 index, qint32 firstNumber, qint32 secondNumber)
{
    clear(false);
    QVector <double> pointsFirst;
    for (qint32 i = 0; i < measuresMap[firstNumber][index].error.size(); ++i)
    {
        pointsFirst.append(i);
    }

    QVector <double> pointsSecond;
    for (qint32 i = 0; i < measuresMap[secondNumber][index].error.size(); ++i)
    {
        pointsSecond.append(i);
    }

    for (auto &i : measuresMap.keys())
    {
        plotter.setPlot(plots[TimeGraph]);
        plotter.setTitle("Times");
        plotter.setGraphName(QString::number(i));
        plotter.setCustomColor(colorMap[i]);
        plotter.addDefaultGraph(measuresMap[i][index].timeDelta, pointsFirst);

        plotter.setPlot(plots[ErrorGraph]);
        plotter.setTitle("Errors");
        plotter.setGraphName(QString::number(i));
        plotter.setCustomColor(colorMap[i]);
        plotter.addDefaultGraph(measuresMap[i][index].error, pointsFirst);
    }
    for (auto& j : plotSyncLines)
    {
        QCustomPlot* p = j->parentPlot();
        j->point1->setCoords(measuresMap[firstNumber][index].time.toMSecsSinceEpoch(), p->yAxis->range().lower);
        j->point2->setCoords(measuresMap[firstNumber][index].time.toMSecsSinceEpoch(), p->yAxis->range().upper);
        p->replot();
    }
}




ApproximationVisualizer::~ApproximationVisualizer()
{
    delete picture;
    delete scenePicture;
}

void ApproximationVisualizer::drawInfoPlaneFullPicture(Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam,
                                                       Calibration::Position2D& XYShift, Calibration::Position& pShift)
{
    drawZones(eOr, cam);
    if (Calibration::GetXYfromXYZ(eOr, cam, pShift, XYShift))
    {
        QPixmap cross("://resources/cross2.png");
        QGraphicsPixmapItem* crossItem =  new QGraphicsPixmapItem(picture);
        crossItem->setPixmap(cross.scaled(cross.width() / 2, cross.height() / 2));
        crossItem->setPos(XYShift.X - cross.width() / 4 , XYShift.Y - cross.height() / 4);
    }
}

QImage ApproximationVisualizer::makeFullPicture(BallApproximator &approx)
{
    scenePicture->removeItem(picture);
    delete picture;
    picture = new QGraphicsPixmapItem(fullPixmap);
    scenePicture->addItem(picture);



    double pos[3], v[3], a[3];
    approx.rotateMovementParameters(pos, v ,a);
    double xp[3];
    approx.getXNonLinearParameters(xp);
    double yp[3];
    approx.getYNonLinearParameters(yp);
    double zp[3];
    approx.getZNonLinearParameters(zp);


    Calibration::Position p, shadowP, pShift;
    approx.calculatePhysicsParameters(tBegin, tEnd, T, vBegin,
                                      vEnd, dxNoRot, dzNoRot, zBegin,
                                      xBegin, rot, W, tFarZone);
    pShift.X = rot[0];
    pShift.Y = rot[1];
    pShift.Z = rot[2];
    QFont font("Ubuntu", 20);
    const qint32 textPosX = 450;
    QGraphicsTextItem * textItemVBegin = new QGraphicsTextItem(picture);
    textItemVBegin->setPlainText(QString::number(vBegin * metersToMiles, 'f', 1) + QString(" mph"));
    textItemVBegin->setPos(textPosX,230);
    textItemVBegin->setFont(font);

    QGraphicsTextItem* textItemVEnd = new QGraphicsTextItem(picture);
    textItemVEnd->setPlainText(QString::number(vEnd * metersToMiles, 'f', 1) + " mph");
    textItemVEnd->setPos(textPosX,280);
    textItemVEnd->setFont(font);

    QGraphicsTextItem * textItemVp = new QGraphicsTextItem(picture);
    textItemVp->setPlainText(QString::number(vBegin * metersToMiles, 'f', 1) + " mph");
    textItemVp->setPos(textPosX, 330);
    textItemVp->setFont(font);

    QGraphicsTextItem * textItemDxNoRot = new QGraphicsTextItem(picture);
    textItemDxNoRot->setPlainText(QString::number(abs(dzNoRot) * 100, 'f', 1) + " cm");
    textItemDxNoRot->setPos(textPosX, 400);
    textItemDxNoRot->setFont(font);

    QGraphicsTextItem * textItemDzNoRot = new QGraphicsTextItem(picture);
    textItemDzNoRot->setPlainText(QString::number(abs(dxNoRot) * 100, 'f', 1) + " cm");
    textItemDzNoRot->setPos(textPosX, textPosX);
    textItemDzNoRot->setFont(font);

    QGraphicsTextItem * textItemTrueSpin = new QGraphicsTextItem(picture);
    textItemTrueSpin->setPlainText(QString::number(W, 'f', 1) + " rpm");
    textItemTrueSpin->setPos(textPosX,500);
    textItemTrueSpin->setFont(font);

    QGraphicsTextItem * textItemZBegin = new QGraphicsTextItem(picture);
    textItemZBegin->setPlainText(QString::number(zBegin, 'f', 1) + " m");
    textItemZBegin->setPos(textPosX, 570);
    textItemZBegin->setFont(font);

    QGraphicsTextItem * textItemXBegin = new QGraphicsTextItem(picture);
    textItemXBegin->setPlainText(QString::number(abs(xBegin), 'f', 1) + " m");
    textItemXBegin->setPos(textPosX,620);
    textItemXBegin->setFont(font);

    QGraphicsTextItem * textItemExtension = new QGraphicsTextItem(picture);
    textItemExtension->setPlainText(QString("1.88") + " m");
    textItemExtension->setPos(textPosX, 670);
    textItemExtension->setFont(font);


    Calibration::ExteriorOr eOr;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cam;
    CalibrationAdjustHelper::readCurrentCalibrationParameters(9999, "calibrate/",
                                                              eOr, cam, false, 1920, 1080);
    Calibration::Position2D XYShadow, XYShift;

    QPixmap ball("://resources/baseball_PNG18979.png");
    QPixmap shadow("://resources/shadow.png");
    double step = T / (ballCount - 1);

    bool drawPlane = false;
    for (qint32 i = 0; i < ballCount; ++i)
    {
        double t2 = pow(tBegin + step * i, 2);
        double t = tBegin + step * i;
        p.X = xp[0] + xp[1] * (t) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (t) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (t) + zp[2] * t2;
        shadowP.X = p.X;
        shadowP.Y = p.Y;
        shadowP.Z = 0;

        double mInit[3][3] {{1, 0 , 0}, {0, 1, 0}, {0, 0, 1}};
        double mRot[3][3];
        double pRot[3];
        double pNoRot[3] {p.X, p.Y, p.Z};
        BOKZMath::rotateOZ(45.0 * BOKZMath::degreesToRad, mInit, mRot);
        multMatrixVector(mRot, pNoRot, pRot);

        Calibration::Position2D dXY;
        if (Calibration::GetXYfromXYZ(eOr, cam, p, dXY))
        {
            Calibration::GetXYfromXYZ(eOr,cam, shadowP, XYShadow);

            qint32 facticalBallSize = ballSize / (pRot[1] + 2) * coeff;

            QGraphicsPixmapItem* sItem = new QGraphicsPixmapItem(picture);
            sItem->setPixmap(shadow.scaled(facticalBallSize , facticalBallSize / 2));
            sItem->setPos(XYShadow.X - facticalBallSize / 2, XYShadow.Y - facticalBallSize / 4);
            if (i == ballCount - 1)
            {
                drawPlane = true;
                drawInfoPlaneFullPicture(eOr, cam, XYShift, pShift);
            }
            QGraphicsPixmapItem* cItem = new QGraphicsPixmapItem(picture);
            cItem->setPixmap(ball.scaled(facticalBallSize , facticalBallSize));
            cItem->setPos(dXY.X - facticalBallSize / 2, dXY.Y -facticalBallSize / 2);
        }
    }

    if (!drawPlane)
    {
        drawInfoPlaneFullPicture(eOr, cam, XYShift, pShift);
    }

    scenePicture->setSceneRect(scenePicture->itemsBoundingRect());
    QImage image = QImage (scenePicture->sceneRect().size().toSize(), QImage::Format_ARGB32);
    image.fill(Qt::transparent);

    QPainter painter(&image);
    scenePicture->render(&painter);
    //image = image.scaled(scenarioWidth, scenarioHeight);
    return image;
}

using namespace StrikeZone;
void ApproximationVisualizer::drawZones(Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam)
{
    double widthInit = 0;

    QPixmap zoneFar("://resources/zone_far.png");

    Calibration::Position2D zonePosLeft2D;
    Calibration::Position zonePosLeft;

    double mInit [3][3] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double mRot[3][3];
    rotateOZ(-45.0 * degreesToRad, mInit, mRot);

    double zone[3], rotZone[3];
    zone[0] = widthInit - width;
    zone[1] = farZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    zonePosLeft.X = rotZone[0];
    zonePosLeft.Y = rotZone[1];
    zonePosLeft.Z = rotZone[2];

    Calibration::GetXYfromXYZ(eOr, cam, zonePosLeft, zonePosLeft2D);
    Calibration::Position2D zonePosRight2D;
    Calibration::Position zonePosRight;

    zone[0] = widthInit + width;
    zone[1] = farZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    zonePosRight.X = rotZone[0];
    zonePosRight.Y = rotZone[1];
    zonePosRight.Z = rotZone[2];

    Calibration::GetXYfromXYZ(eOr, cam, zonePosRight, zonePosRight2D);

    double xZoneFar = zonePosLeft2D.X;
    double yZoneFar = zonePosLeft2D.Y;
    qint32 scaleX = zonePosRight2D.X - zonePosLeft2D.X;
    qint32 scaleY = zonePosRight2D.Y - zonePosLeft2D.Y;
    zoneFar = zoneFar.scaled(scaleX, scaleY);
    QGraphicsPixmapItem* zoneFarItem =  new QGraphicsPixmapItem(picture);
    zoneFarItem->setPixmap(zoneFar);
    zoneFarItem->setPos(xZoneFar, yZoneFar);


    QPixmap fone("://resources/fone.png");
    zone[0] = widthInit - width;
    zone[1] = closeZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    zonePosLeft.X = rotZone[0];
    zonePosLeft.Y = rotZone[1];
    zonePosLeft.Z = rotZone[2];
    Calibration::GetXYfromXYZ(eOr, cam, zonePosLeft, zonePosLeft2D);

    zone[0] = widthInit + width;
    zone[1] = closeZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    zonePosRight.X = rotZone[0];
    zonePosRight.Y = rotZone[1];
    zonePosRight.Z = rotZone[2];
    Calibration::GetXYfromXYZ(eOr, cam, zonePosRight, zonePosRight2D);


    double xFone = zonePosLeft2D.X;
    double yFone = zonePosLeft2D.Y;
    scaleX = zonePosRight2D.X - zonePosLeft2D.X;
    scaleY = zonePosRight2D.Y - zonePosLeft2D.Y;
    fone = fone.scaled(scaleX, scaleY);
    QGraphicsPixmapItem* foneItem =  new QGraphicsPixmapItem(picture);
    foneItem->setPixmap(fone);
    foneItem->setPos(xFone, yFone);

    QPixmap zoneClose("://resources//zone_close.png");
    double xZoneClose = zonePosLeft2D.X;
    double yZoneClose = zonePosLeft2D.Y;
    zoneClose = zoneClose.scaled(scaleX, scaleY);
    QGraphicsPixmapItem* zoneCloseItem =  new QGraphicsPixmapItem(picture);
    zoneCloseItem->setPixmap(zoneClose);
    zoneCloseItem->setPos(xZoneClose, yZoneClose);
}

void ApproximationVisualizer::drawLastShortPicture(bool& strike, double& vBegin, Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam,
                                                   Calibration::Position2D& XYShift, Calibration::Position& pShift)
{
    QPixmap resFone("://resources/zone_res.png");
    auto resFoneScaled = resFone.scaled(resFone.width() , resFone.height() /*+ 120*/);
    double xFoneRes = 156;
    double yFoneRes = 188;
    QGraphicsPixmapItem* resFoneItem =  new QGraphicsPixmapItem(picture);
    resFoneItem->setPixmap(resFoneScaled);
    resFoneItem->setPos(xFoneRes, yFoneRes);

    QPixmap pitch("://resources/pitching.png");
    double xTitle = 215;
    double yTitle = 220;
    QGraphicsPixmapItem* pitchItem =  new QGraphicsPixmapItem(picture);
    pitchItem->setPixmap(pitch);
    pitchItem->setPos(xTitle, yTitle);

    double xRes;
    double yRes;
    QGraphicsPixmapItem* resItem = new QGraphicsPixmapItem(picture);
    QPixmap res;
    if (strike)
    {
        xRes = 210;
        yRes = 350;
        res = QPixmap("://resources/strike_res.png");
    }
    else
    {
        xRes = 280;
        yRes = 350;
        res = QPixmap("://resources/ball_res.png");
    }

    resItem->setPixmap(res);
    resItem->setPos(xRes, yRes);

    double xSpeed = 180;
    double ySpeed = 435;
    QGraphicsTextItem * textItemVBegin = new QGraphicsTextItem(picture);
    textItemVBegin->setDefaultTextColor(QColor(255, 255, 255));
    textItemVBegin->setPlainText(QString::number(vBegin * metersToMiles, 'f', 1));
    textItemVBegin->setPos(xSpeed, ySpeed);
    auto effect = new QGraphicsDropShadowEffect(textItemVBegin);
    effect->setColor(QColor(0, 0, 0));
    textItemVBegin->setGraphicsEffect(effect);

    QFont font("Segoe UI", 100);
    font.setBold(true);
    font.setItalic(true);

    textItemVBegin->setFont(font);

    QPixmap speed("://resources/speed.png");
    double xSpeedTitle = 480;
    double ySpeedTitle = 510;
    QGraphicsPixmapItem* speedTitleItem =  new QGraphicsPixmapItem(picture);
    speedTitleItem->setPixmap(speed);
    speedTitleItem->setPos(xSpeedTitle, ySpeedTitle);


    drawZones(eOr, cam);
    if (Calibration::GetXYfromXYZ(eOr, cam, pShift, XYShift))
    {
        QPixmap cross("://resources/cross2.png");
        QGraphicsPixmapItem* crossItem =  new QGraphicsPixmapItem(picture);
        crossItem->setPixmap(cross.scaled(cross.width() / 2, cross.height() / 2));
        crossItem->setPos(XYShift.X - cross.width() / 4 , XYShift.Y - cross.height() / 4);
    }
}

QImage ApproximationVisualizer::makeShortPicture(BallApproximator &approx, QString info)
{
    scenePicture->removeItem(picture);
    delete picture;
    picture = new QGraphicsPixmapItem(shortPixmap);
    scenePicture->addItem(picture);

    double pos[3], v[3], a[3];
    approx.rotateMovementParameters(pos, v ,a);
    double xp[3];
    approx.getXNonLinearParameters(xp);
    double yp[3];
    approx.getYNonLinearParameters(yp);
    double zp[3];
    approx.getZNonLinearParameters(zp);

    Calibration::Position p, shadowP, pShift;
    bool strike = approx.calculatePhysicsParameters(tBegin, tEnd, T, vBegin,
                                                    vEnd, dxNoRot, dzNoRot, zBegin,
                                                    xBegin, rot, W, tFarZone);
    pShift.X = rot[0];
    pShift.Y = rot[1];
    pShift.Z = rot[2];

    Calibration::ExteriorOr eOr;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cam;
    CalibrationAdjustHelper::readCurrentCalibrationParameters(9999, "calibrate/",
                                                              eOr, cam, false, 1920, 1080);
    Calibration::Position2D XYShadow, XYShift;


    double xSpeed = 900;
    double ySpeed = 435;
    QGraphicsTextItem * textItemVBegin = new QGraphicsTextItem(picture);
    textItemVBegin->setDefaultTextColor(QColor(255, 255, 255));
    textItemVBegin->setPlainText(info);
    textItemVBegin->setPos(xSpeed, ySpeed);
    QFont font("Segoe UI", 50);
    font.setBold(true);
    font.setItalic(true);

    textItemVBegin->setFont(font);

    QPixmap ball("://resources/baseball_PNG18979.png");
    QPixmap shadow("://resources/shadow.png");
    double step = T / (ballCount - 1);
    bool drawInfoPlane = false;
    for (qint32 i = 0; i < ballCount; ++i)
    {
        double t2 = pow(tBegin + step * i, 2);
        double t = tBegin + step * i;
        p.X = xp[0] + xp[1] * (t) + xp[2] * t2;
        p.Y = yp[0] + yp[1] * (t) + yp[2] * t2;
        p.Z = zp[0] + zp[1] * (t) + zp[2] * t2;
        shadowP.X = p.X;
        shadowP.Y = p.Y;
        shadowP.Z = 0;

        double mInit[3][3] {{1, 0 , 0}, {0, 1, 0}, {0, 0, 1}};
        double mRot[3][3];
        double pRot[3];
        double pNoRot[3] {p.X, p.Y, p.Z};
        BOKZMath::rotateOZ(45.0 * BOKZMath::degreesToRad, mInit, mRot);
        multMatrixVector(mRot, pNoRot, pRot);

        Calibration::Position2D dXY;
        if (Calibration::GetXYfromXYZ(eOr, cam, p, dXY))
        {
            Calibration::GetXYfromXYZ(eOr,cam, shadowP, XYShadow);

            qint32 facticalBallSize = ballSize / (pRot[1] + 2) * coeff;

            QGraphicsPixmapItem* sItem = new QGraphicsPixmapItem(picture);
            sItem->setPixmap(shadow.scaled(facticalBallSize , facticalBallSize / 2));
            sItem->setPos(XYShadow.X - facticalBallSize / 2, XYShadow.Y - facticalBallSize / 4);
            if (i == ballCount - 1)
            {
                drawInfoPlane = true;
                drawLastShortPicture(strike, vBegin, eOr, cam, XYShift, pShift);
            }
            QGraphicsPixmapItem* cItem = new QGraphicsPixmapItem(picture);
            cItem->setPixmap(ball.scaled(facticalBallSize , facticalBallSize));
            cItem->setPos(dXY.X - facticalBallSize / 2, dXY.Y -facticalBallSize / 2);
        }
    }
    if (!drawInfoPlane)
    {
        drawLastShortPicture(strike, vBegin, eOr, cam, XYShift, pShift);
    }

    scenePicture->setSceneRect(scenePicture->itemsBoundingRect());
    QImage image = QImage (scenePicture->sceneRect().size().toSize(), QImage::Format_ARGB32);
    image.fill(Qt::transparent);

    QPainter painter(&image);
    scenePicture->render(&painter);
    QImage scaled = image.scaled(scenarioWidth, scenarioHeight);
    return scaled;
}

QImage ApproximationVisualizer::addHitInfo(Mat img, HitParameters& hParams)
{
    scenePicture->removeItem(picture);
    delete picture;
    QPixmap px = QPixmap::fromImage(QImage((unsigned char*) img.data, img.cols, img.rows, QImage::Format_RGB888));
    picture = new QGraphicsPixmapItem(px);
    scenePicture->addItem(picture);


    QPixmap resFone("://resources/zone_res.png");
    auto resFoneScaled = resFone.scaled(resFone.width() / 1.5, resFone.height() / 1.5 /*+ 80*/);
    double xFoneRes = 785;
    double yFoneRes = 125;
    QGraphicsPixmapItem* resFoneItem =  new QGraphicsPixmapItem(picture);
    resFoneItem->setPixmap(resFoneScaled);
    resFoneItem->setPos(xFoneRes, yFoneRes);

    QPixmap batt("://resources/batting.png");
    auto battScaled = batt.scaled(batt.width() / 1.5, batt.height() / 1.5);
    double xTitle = 839;
    double yTitle = 145;
    QGraphicsPixmapItem* pitchItem =  new QGraphicsPixmapItem(picture);
    pitchItem->setPixmap(battScaled);
    pitchItem->setPos(xTitle, yTitle);

    QFont font("Segoe UI", 66);
    font.setBold(true);
    font.setItalic(true);
    QGraphicsTextItem * textItemVBegin = new QGraphicsTextItem(picture);
    textItemVBegin->setDefaultTextColor(QColor(255, 255, 255));
    textItemVBegin->setPlainText(QString::number(hParams.initSpeed, 'f', 0));
    double xSpeed = 809;
    double ySpeed = 195;
    textItemVBegin->setPos(xSpeed, ySpeed);
    auto effect = new QGraphicsDropShadowEffect(textItemVBegin);
    effect->setColor(QColor(0, 0, 0));
    textItemVBegin->setGraphicsEffect(effect);
    textItemVBegin->setFont(font);

    QPixmap speed("://resources/speed.png");
    auto speedScaled = speed.scaled(speed.width() / 1.5, speed.height() / 1.5);
    double xSpeedTitle = 980;
    double ySpeedTitle = 245;
    QGraphicsPixmapItem* speedTitleItem =  new QGraphicsPixmapItem(picture);
    speedTitleItem->setPixmap(speedScaled);
    speedTitleItem->setPos(xSpeedTitle, ySpeedTitle);

    QGraphicsTextItem * textItemAngle = new QGraphicsTextItem(picture);
    textItemAngle->setDefaultTextColor(QColor(255, 255, 255));
    textItemAngle->setPlainText(QString::number(hParams.angle, 'f', 0) + "°");
    double xAngle = 800;
    double yAngle = 290;
    textItemAngle->setPos(xAngle, yAngle);
    effect = new QGraphicsDropShadowEffect(textItemAngle);
    effect->setColor(QColor(0, 0, 0));
    textItemAngle->setGraphicsEffect(effect);
    textItemAngle->setFont(font);

    //    QGraphicsTextItem * textItemDistance = new QGraphicsTextItem(picture);
    //    textItemDistance->setDefaultTextColor(QColor(255, 255, 255));
    //    textItemDistance->setPlainText(QString::number(hParams.distance, 'f', 0));
    //    double xDistance = 800;
    //    double yDistance = 385;
    //    textItemDistance->setPos(xDistance, yDistance);
    //    effect = new QGraphicsDropShadowEffect(textItemDistance);
    //    effect->setColor(QColor(0, 0, 0));
    //    textItemDistance->setGraphicsEffect(effect);
    //    textItemDistance->setFont(font);

    //    QPixmap meters("://resources/meters.png");
    //    auto metersScaled = meters.scaled(meters.width() / 1.5, meters.height() / 1.5);
    //    double xMetersTitle = 980;
    //    double yMetersTitle = 450;
    //    QGraphicsPixmapItem* metersTitleItem =  new QGraphicsPixmapItem(picture);
    //    metersTitleItem->setPixmap(metersScaled);
    //    metersTitleItem->setPos(xMetersTitle, yMetersTitle);

    scenePicture->setSceneRect(scenePicture->itemsBoundingRect());
    QImage image = QImage (scenePicture->sceneRect().size().toSize(), QImage::Format_ARGB32);
    image.fill(Qt::transparent);

    QPainter painter(&image);
    scenePicture->render(&painter);
    QImage scaled = image.scaled(scenarioWidth, scenarioHeight);
    return scaled;


}

void ApproximationVisualizer::readDrawTracerDebugData(const QString& fVideo, Calibration::ExteriorOr& EOFirstCamera,
                                                      Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                                                      QVector <Calibration::Position>& firstVecs, QVector <double>& firstTime, QVector <double>& videoTimes)
{
    qint32 num;
    if (fVideo.contains("3850")) num = 3850;
    else num = 4510;
    QString timesPath = fVideo;
    timesPath.remove(".avi").append(".txt");
    QFile fTimePath(timesPath);
    QVector <quint64> fTimes;
    if (fTimePath.open(QIODevice::ReadOnly))
    {
        QTextStream in (&fTimePath);
        QString line;
        while (in.readLineInto(&line))
        {
            videoTimes.append((double)line.toLongLong() / 10000000.0);
        }
    }
    QString fCoordPath = timesPath;
    fCoordPath.remove(".txt").append("_coord.txt");
    QFile fCoordFile(fCoordPath);
    QVector <QStringList> points;
    Calibration::RayAndPoint rp;
    Calibration::Position2D XYpix_left;

    if (fCoordFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&fCoordFile);
        QString line;

        while (in.readLineInto(&line))
        {
            points.append(line.split("\t"));
        }
    }

    CalibrationAdjustHelper::readCurrentCalibrationParameters(num, "calibrate", EOFirstCamera, cameraFirst, true);
    QVector <Point2f> balls;
    for (qint32 i = 0; i < points.size(); ++i)
    {
        XYpix_left.X = points[i][0].toDouble();
        XYpix_left.Y = points[i][1].toDouble();
        balls.append(Point2f(XYpix_left.X, XYpix_left.Y));
        Calibration::GetRayAndPoint(EOFirstCamera, cameraFirst, XYpix_left, rp);
        firstVecs.append(rp.Vect);
        // qDebug() << rp.Vect.X << rp.Vect.Y  << rp.Vect.Z << QString::number(points[i][2].toDouble()/ divideTime, 'g', 10);
        double time = points[i][2].toDouble();
        if (qFuzzyCompare(points[i][3].toDouble(), 1))
        {
            firstTime.append(time);
        }
        else
        {
            firstTime.append(-time);
        }
    }

    //    auto p1 = balls[balls.size() / 2];
    //    auto p2 = balls[(balls.size() / 2) + 1];
    //    auto p3 = p2 - p1;
    //    auto l  = norm(p3);
    //    auto p4 = p3 / l;
    //    qDebug() << num << p4.x << p4.y;
    //    auto p5 = Point2f(779.591,	408.669);
    //    auto p6 = Point2f(766.168,	400.631);
    //    auto p7 = p6 - p5;
    //    l = norm(p7);
    //    auto p8 = p7 / l;
    //    qDebug() << acos(p4.dot(p8)) * radToDegrees;
}

void ApproximationVisualizer::drawStrikeZoneParallelepiped(Calibration::ExteriorOr& EOFirstCamera,
                                                           Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                                                           QGraphicsItem* parent)
{
    QPolygonF fPlane;
    Calibration::Position p;
    double mInit [3][3] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double mRot[3][3];
    rotateOZ(-45.0 * degreesToRad, mInit, mRot);
    double zone[3], rotZone[3];
    zone[0] = -width;
    zone[1] = farZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p1;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p1);


    zone[0] = -width;
    zone[1] = closeZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p2;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2);

    zone[0] = -width;
    zone[1] = closeZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p3;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p3);


    zone[0] = -width;
    zone[1] = farZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p4;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p4);
    fPlane << QPointF(p1.X, p1.Y) << QPointF(p2.X, p2.Y) << QPointF(p3.X, p3.Y) << QPointF(p4.X, p4.Y) << QPointF(p1.X, p1.Y);


    QPolygonF sPlane;
    zone[0] = width;
    zone[1] = farZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p5;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p5);

    zone[0] = width;
    zone[1] = farZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p6;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p6);
    sPlane << QPointF(p1.X, p1.Y) << QPointF(p4.X, p4.Y) << QPointF(p5.X, p5.Y) << QPointF(p6.X, p6.Y) << QPointF(p1.X, p1.Y);


    QPolygonF tPlane;
    zone[0] = width;
    zone[1] = closeZoneY;
    zone[2] = maxHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p7;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p7);
    zone[0] = width;
    zone[1] = closeZoneY;
    zone[2] = minHeight;
    multMatrixVector(mRot, zone, rotZone);
    p.X = rotZone[0];
    p.Y = rotZone[1];
    p.Z = rotZone[2];
    Calibration::Position2D p8;
    Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p8);
    tPlane << QPointF(p5.X, p5.Y) << QPointF(p6.X, p6.Y) << QPointF(p7.X, p7.Y) << QPointF(p8.X, p8.Y) << QPointF(p5.X, p5.Y);


    QPolygonF frPlane;
    frPlane << QPointF(p7.X, p7.Y) << QPointF(p8.X, p8.Y) << QPointF(p3.X, p3.Y) << QPointF(p2.X, p2.Y) << QPointF(p7.X, p7.Y);

    QPen pen = QPen(QBrush(QColor(255, 0, 0, 200)), 2, Qt::SolidLine, Qt::SquareCap);
    QGraphicsPolygonItem* fPoly = new QGraphicsPolygonItem(parent);
    fPoly->setPolygon(fPlane);
    fPoly->setPen(pen);
    QGraphicsPolygonItem* sPoly = new QGraphicsPolygonItem(parent);
    sPoly->setPolygon(sPlane);
    sPoly->setPen(pen);
    QGraphicsPolygonItem* tPoly = new QGraphicsPolygonItem(parent);
    tPoly->setPolygon(tPlane);
    tPoly->setPen(pen);
    QGraphicsPolygonItem* frPoly = new QGraphicsPolygonItem(parent);
    frPoly->setPolygon(frPlane);
    frPoly->setPen(pen);
}

void ApproximationVisualizer::runRepeatStreamInternal()
{
    repeatThreadRun = 1;
    if (!outWindow.isOpened())
    {
        Q_ASSERT(outWindow.open("appsrc ! videoconvert ! d3dvideosink sync=false async=false", CV_CAP_GSTREAMER, 30, Size(scenarioWidth, scenarioHeight)));
    }
    RepeatVisualizeData repeatData;
    bool newRepeat = false;
    qDebug() << "SCENARIO STARTED";
    qint32 handledFrames = 0;
    RepeatVisualizeState prevState;
    while (repeatThreadRun && outWindow.isOpened())
    {

        switch (currentRepeatState)
        {
        case ShowInit:
        {
            handledFrames = 0;
            newRepeat = false;
            QMutexLocker locker(&scenarioMutex);
            outWindow.write(scenarioData.initPicture);
            locker.unlock();
            QThread::msleep(1000. / 30.);
            break;
        }
        case ShowResultPicture:
        {
            handledFrames = 0;
            newRepeat = false;
            QMutexLocker locker(&scenarioMutex);
            outWindow.write(scenarioData.resultPicture);
            locker.unlock();
            QThread::msleep(1000. / 30.);
            break;
        }
        case ShowRepeat:
        {
            if (!newRepeat)
            {
                qDebug() << "NEW REPEAT" << scenarioData.repeats.size();
                repeatData = RepeatVisualizeData();
                newRepeat = true;
            }
            QMutexLocker locker(&scenarioMutex);
            if (scenarioData.repeats.size() > 1)
            {
                double corrFirst = scenarioData.repeats.first().first.corr;
                double corrSecond = scenarioData.repeats.last().first.corr;
                qDebug() << "CHECK CONFLICT" << corrFirst << corrSecond;
                if (qFuzzyCompare(corrFirst, -1) || qFuzzyCompare(corrSecond, -1))
                {
                    continue;
                }
                else
                {
                    double timeFirst = scenarioData.repeats.first().first.startTime;
                    double timeSecond = scenarioData.repeats.last().first.startTime;
                    if (scenarioData.approx.isNull()
                            || abs(timeFirst - timeSecond) > 1)
                    {
                        qDebug() << "INCORRECT REPEATS" << timeFirst << timeSecond;
                        scenarioData.repeats.clear();
                        currentRepeatState = ShowInit;
                        continue;
                    }
                    qDebug() << "RESOLVE CONFLICT" << corrFirst << corrSecond << scenarioData.repeats.keys().first() << scenarioData.repeats.keys().last();
                    if (corrFirst > corrSecond)
                    {
                        scenarioData.repeats.remove(scenarioData.repeats.keys().last());
                    }
                    else
                    {
                        scenarioData.repeats.remove(scenarioData.repeats.keys().first());
                    }
                    scenarioData.corrConflictResolved = true;
                    setRepeatCameraNumber(scenarioData.repeats.keys().first());
                }
            }

            if (!scenarioData.repeats.isEmpty())
            {
                auto& v = scenarioData.repeats.first().second;
                if (!v.isEmpty() && scenarioData.corrConflictResolved)
                {
                    auto& v = scenarioData.repeats.first().second;
                    QPair <cv::Mat, double> pair = v.first();
                    Mat raw;
                    pair.first.copyTo(raw); // вот тут сделать просто raw = pair.first
                    v.removeFirst();
                    cvtColor(raw, raw, CV_BGR2RGB);
                    QScopedPointer <QElapsedTimer> t (new QElapsedTimer());
                    t->start();
                    QPixmap px = QPixmap::fromImage(QImage((unsigned char*) raw.data, raw.cols, raw.rows, QImage::Format_RGB888));
                    Mat m = drawBallTracer(repeatData, px, scenarioData.EOcamera, scenarioData.camera,
                                           *scenarioData.approx, pair.second,
                                           scenarioData.repeats.first().first.initTime);
                    Mat preScaled, scaled;
                    m(scenarioData.scaleRects[scenarioData.repeats.keys().first()]).copyTo(preScaled);
                    ++handledFrames;
                    if (handledFrames == scenarioData.repeats.first().first.frameCount)
                    {
                        currentRepeatState = ShowResultPicture;
                    }
                    locker.unlock();
                    cv::resize(preScaled, scaled, Size(scenarioWidth, scenarioHeight));
                    qint64 elapsed = t->elapsed();
                    outWindow.write(scaled);
                    qint32 needToWait = (1000. / 30.) - elapsed;

                    if (needToWait > 0)
                    {
                        QThread::msleep(needToWait);
                    }


                }
            }
            break;
        }
        }

        prevState = currentRepeatState;
    }
    if (repeatThreadRun && ! outWindow.isOpened())
    {
        qDebug() << "TROUBLE in scenario thread";
    }
}

cv::Mat ApproximationVisualizer::drawBallTracer(RepeatVisualizeData& repeatData, QPixmap& px, Calibration::ExteriorOr& EOFirstCamera,
                                                Calibration::SpacecraftPlatform::CAMERA::CameraParams& cameraFirst,
                                                BallApproximator& approx, double videoTime, double initTime)
{
    QGraphicsPixmapItem* repeatPicture = new QGraphicsPixmapItem(px);
    QScopedPointer <QGraphicsScene> sceneRepeat(new QGraphicsScene());
    sceneRepeat->addItem(repeatPicture);
    QGraphicsPathItem* pathItem = new QGraphicsPathItem(repeatPicture);
    //QGraphicsPathItem* areaPathItem = new QGraphicsPathItem(repeatPicture);
    //QGraphicsPathItem* linesPathItem = new QGraphicsPathItem(repeatPicture);
    QGraphicsPathItem* bottomPathItem = new QGraphicsPathItem(repeatPicture);

    double point[3];
    double delta = 166200.0 / 10000000.0;
    Calibration::Position2D p2d;
    Calibration::Position p;

    drawStrikeZoneParallelepiped(EOFirstCamera, cameraFirst, repeatPicture);
    //qDebug() << "2";
    double startTime = approx.getTIN();
    videoTime = videoTime - startTime;
    initTime = abs(initTime) - startTime;
    if (abs(videoTime - initTime) < delta && !repeatData.startPointFound)
    {
        repeatData.startPointFound = true;
        approx.getPointAt(videoTime, point);
        p.X = point[0];
        p.Y = point[1];
        p.Z = point[2];
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
        repeatData.path.moveTo(p2d.X, p2d.Y);

        p.X = point[0];
        p.Y = point[1];
        p.Z = 0;
        repeatData.p2dPrevUp = p2d;
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
        repeatData.p2dPrevDown = p2d;
    }

    if (repeatData.startPointFound && !repeatData.isLastPoint)
    {
        double tFactical = videoTime;

        approx.getPointAt(tFactical, point);

        double mInit [3][3] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        double mRot[3][3];
        rotateOZ(45.0 * degreesToRad, mInit, mRot);
        double rot[3];
        multMatrixVector(mRot, point, rot);

        if (rot[1] < closeZoneY)
        {
            approx.getPointAt(tEnd, point);
            p.X = point[0];
            p.Y = point[1];
            p.Z = point[2];
            Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
            repeatData.closeZoneEllipse.setX(p2d.X - 5);
            repeatData.closeZoneEllipse.setY(p2d.Y - 5);
            repeatData.closeZoneEllipse.setWidth(10);
            repeatData.closeZoneEllipse.setHeight(10);
        }
        if (rot[1] < farZoneY)
        {
            approx.getPointAt(tFarZone, point);
            p.X = point[0];
            p.Y = point[1];
            p.Z = point[2];
            Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
            repeatData.farZoneEllipse.setX(p2d.X - 5);
            repeatData.farZoneEllipse.setY(p2d.Y - 5);
            repeatData.farZoneEllipse.setWidth(10);
            repeatData.farZoneEllipse.setHeight(10);
        }

        if (rot[1] < closeZoneY)
        {
            repeatData.isLastPoint = true;
            tFactical = tEnd;
            approx.getPointAt(tFactical, point);
        }

        p.X = point[0];
        p.Y = point[1];
        p.Z = point[2];
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
        repeatData.path.lineTo(p2d.X, p2d.Y);
        p.X = point[0];
        p.Y = point[1];
        p.Z = 0;
        QPolygonF pl;
        auto plUp = QPointF(repeatData.p2dPrevUp.X, repeatData.p2dPrevUp.Y);
        auto plDown = QPointF(repeatData.p2dPrevDown.X, repeatData.p2dPrevDown.Y);
        repeatData.p2dPrevUp = p2d;
        Calibration::GetXYfromXYZ(EOFirstCamera, cameraFirst, p, p2d);
        repeatData.p2dPrevDown = p2d;
        pl << plUp <<  plDown << QPointF(repeatData.p2dPrevDown.X, repeatData.p2dPrevDown.Y)
           << QPointF(repeatData.p2dPrevUp.X, repeatData.p2dPrevUp.Y);
        repeatData.polygons.append(pl);
        repeatData.areaPath.addPolygon(pl);
        repeatData.areaPath.closeSubpath();
        repeatData.linesPath.moveTo(repeatData.p2dPrevUp.X, repeatData.p2dPrevUp.Y + 10);
        repeatData.linesPath.lineTo(repeatData.p2dPrevDown.X, repeatData.p2dPrevDown.Y);
        repeatData.bottomPath.moveTo(plDown);
        repeatData.bottomPath.lineTo(repeatData.p2dPrevDown.X, repeatData.p2dPrevDown.Y);

    }

    int cRed = 255;
    int cBlue = 206;
    for (qint32 j = 0; j < repeatData.polygons.size(); ++j)
    {
        QGraphicsPolygonItem* poly = new QGraphicsPolygonItem(repeatPicture);
        poly->setPolygon(repeatData.polygons[j]);

        QColor c;
        if (j <= approx.getFirstTime().size() / 2)
        {
            c = QColor(cRed, cBlue, 250, 50);
            cRed -= 5;
            cBlue -= 3;
        }
        else
        {
            if (cRed > UCHAR_MAX)
            {
                cRed = UCHAR_MAX;
            }
            if (cBlue > UCHAR_MAX)
            {
                cBlue = UCHAR_MAX;
            }
            c = QColor(cRed, cBlue, 250, 50);
            cRed += 5;
            cBlue += 3;
        }
        poly->setBrush(QBrush(c));
        poly->setPen(QPen(QBrush(c), 1, Qt::SolidLine, Qt::SquareCap));
    }
    //        areaPathItem->setPen(QPen(QBrush(QColor(135, 206, 250, 20)), 1, Qt::SolidLine, Qt::RoundCap));
    //        areaPathItem->setBrush(QBrush(QColor(135, 206, 250, 50)));
    //        areaPathItem->setPath(areaPath);

    //linesPathItem->setPen(QPen(QBrush(/*grad*/QColor(135, 206, 250, 100)), 20, Qt::SolidLine, Qt::SquareCap));
    //linesPathItem->setPath(linesPath);
    //auto blur = new QGraphicsBlurEffect();
    //blur->setBlurRadius(3);


    bottomPathItem->setPen(QPen(QBrush(QColor(30, 144, 255, 150)), 7, Qt::SolidLine, Qt::RoundCap));
    bottomPathItem->setPath(repeatData.bottomPath);

    pathItem->setPen(QPen(QBrush(QColor(135, 206, 250, 100)), 7, Qt::SolidLine, Qt::RoundCap));
    pathItem->setPath(repeatData.path);


    sceneRepeat->setSceneRect(sceneRepeat->itemsBoundingRect());
    QImage image = QImage (sceneRepeat->sceneRect().size().toSize(), QImage::Format_ARGB32);
    //image.fill(Qt::transparent);

    if (!repeatData.farZoneEllipse.isNull())
    {
        QGraphicsEllipseItem* farZoneItem = new QGraphicsEllipseItem(repeatPicture);
        farZoneItem->setBrush(QBrush(QColor(255, 255, 255, 150)));
        farZoneItem->setPen(QPen(QBrush(QColor(255, 255, 255, 150)), 1, Qt::SolidLine, Qt::SquareCap));
        farZoneItem->setRect(repeatData.farZoneEllipse);
    }
    if (!repeatData.closeZoneEllipse.isNull())
    {
        QGraphicsEllipseItem* closeZoneItem = new QGraphicsEllipseItem(repeatPicture);
        closeZoneItem->setBrush(QBrush(QColor(255, 255, 255, 150)));
        closeZoneItem->setPen(QPen(QBrush(QColor(255, 255, 255, 150)), 1, Qt::SolidLine, Qt::SquareCap));
        closeZoneItem->setRect(repeatData.closeZoneEllipse);
    }

    QPainter painter(&image);
    sceneRepeat->render(&painter);
    Mat m(image.height(), image.width(), CV_8UC4, (void*)image.constBits());
    cvtColor(m, m, CV_RGBA2RGB);
    return m;
}

void ApproximationVisualizer::drawTracerDebug(const QString &fVideo, const QString &sVideo)
{


    Calibration::ExteriorOr EOFirstCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraFirst;
    QVector <Calibration::Position> firstVecs;
    QVector <double> firstTime,  firstVideoTime;
    readDrawTracerDebugData(fVideo, EOFirstCamera, cameraFirst, firstVecs, firstTime, firstVideoTime);

    Calibration::ExteriorOr EOSecondCamera;
    Calibration::SpacecraftPlatform::CAMERA::CameraParams cameraSecond;
    QVector <Calibration::Position> secondVecs;
    QVector <double> secondTime,  secondVideoTime;
    readDrawTracerDebugData(sVideo, EOSecondCamera, cameraSecond, secondVecs, secondTime, secondVideoTime);

    //    QString timesPath = fVideo;
    //    timesPath.remove(".avi").append(".txt");

    QSharedPointer <BallApproximator> approx(new BallApproximator);
    approx->readData(firstVecs, secondVecs, firstTime, secondTime, EOFirstCamera.Point, EOSecondCamera.Point);
    approx->calculateApproximation(QString("games_%1/results/result%2.txt")
                                   .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                   .arg(QTime::currentTime().toString("hh_mm_ss")), true);

    auto v1 = approx->getFirstErrors();
    //qDebug() << "error1" << v1;
    double coef = 0.05;
    bool repeat = false;
    while (v1.size() > 0 && v1.first() > coef)
    {
        firstTime.removeFirst();
        firstVecs.removeFirst();
        v1.removeFirst();
        repeat = true;
    }
    while (v1.size() > 0  && v1.last() > coef)
    {
        firstTime.removeLast();
        firstVecs.removeLast();
        v1.removeLast();
        repeat = true;
    }
    auto v2 = approx->getSecondErrors();
    //qDebug() << "error2" << v2;
    while ( v2.size() > 0  &&  v2.first() > coef)
    {
        secondTime.removeFirst();
        secondVecs.removeFirst();
        v2.removeFirst();
        repeat = true;
    }
    while ( v2.size() > 0 && v2.last() > coef)
    {
        secondTime.removeLast();
        secondVecs.removeLast();
        v2.removeLast();
        repeat = true;
    }
    if (repeat)
    {
        approx->readData(firstVecs, secondVecs, firstTime, secondTime, EOFirstCamera.Point, EOSecondCamera.Point);
        approx->calculateApproximation(QString("games_%1/results/result%2.txt")
                                       .arg(QDate::currentDate().toString("dd_MM_yyyy"))
                                       .arg(QTime::currentTime().toString("hh_mm_ss")), true);

    }


    cv::VideoCapture capFirst = cv::VideoCapture(fVideo.toStdString());
    qint32 length = capFirst.get(cv::CAP_PROP_FRAME_COUNT);
    double initTime = firstTime.first();
    //    bool startPointFound = false;
    //    double point[3];
    //    double delta = 166200.0 / 10000000.0;
    //Calibration::Position2D p2d, p2dPrevUp, p2dPrevDown;
    //Calibration::Position p;
    QPainterPath path;
    QPainterPath bottomPath;
    QPainterPath areaPath;
    QPainterPath linesPath;
    QVector <QPolygonF> polygons;

    //cv::VideoWriter outWindow("appsrc ! videoconvert ! d3dvideosink sync=false async=false", CV_CAP_GSTREAMER, 30, Size(1920, 1080));
    outWindow.open("appsrc ! videoconvert ! d3dvideosink sync=false async=false", CV_CAP_GSTREAMER, 30, Size(scenarioWidth, scenarioHeight));
    //cv::VideoWriter outWindow("appsrc ! videoconvert ! glimagesink sync=false async=false", CV_CAP_GSTREAMER, 30, Size(1920, 1080));
    if (outWindow.isOpened())
    {
        qDebug() << "QQ";
    }

    cv::VideoWriter video;
    video.open(QString("games_%1/test_video.avi")
               .arg(QDate::currentDate().toString("dd_MM_yyyy")).toStdString(), CV_FOURCC('X','V','I','D'), 60, cv::Size(1920, 1080));
    //bool isLastPoint = false;
    QVector <cv::Mat> data;
    runRepeatStream();
    QThread::msleep(1000);
    setCurrentRepeatState(ShowRepeat);
    setRepeatCameraNumber(3850);
    setRepeatInitTime(3850, initTime);
    setLastApproximation(approx);
    for (qint32 i = 0; i < length; ++i)
    {
        cv::Mat frame;
        capFirst >> frame;
        appendRepeatFrame(3850, frame, firstVideoTime[i]);
        //data.append(frame);
        // qDebug() << "qq1";
        //QThread::msleep(100);


        qDebug() << "qq2";
        //        qDebug() << calculateBatterPositionCorr(3850, frame);
        //        //break;
        //        QElapsedTimer t;
        //        t.start();
        //        cvtColor(frame, frame, CV_BGR2RGB);
        //        QPixmap px = QPixmap::fromImage(QImage((unsigned char*) frame.data, frame.cols, frame.rows, QImage::Format_RGB888));
        //        auto m = drawBallTracer(px, EOFirstCamera, cameraFirst, approx, firstVideoTime[i], initTime);
        //        qDebug() << t.elapsed();
        //video.write(m);
    }
    //    qint32 j = 1;
    //    while (true)
    //    {

    //        for (auto& i : data)
    //        {
    //            QThread::msleep(180u / j);

    //            outWindow.write(i);

    //        }
    //        j += 1;
    //        if (j == 6) break;
    //    }

    //    // outWindow.release();
    //    cv::VideoCapture capSecond = cv::VideoCapture(sVideo.toStdString());
    //    length = capSecond.get(cv::CAP_PROP_FRAME_COUNT);
    //    cv::Mat frame;
    //    capSecond >> frame;
    //    qDebug() << calculateBatterPositionCorr(4510, frame);
}

void ApproximationVisualizer::clearCalibGraphs()
{
    if (!plots.isEmpty())
    {
        plots[CalibFirstGraph]->clearGraphs();
        plots[CalibSecondGraph]->clearGraphs();
    }
}

double ApproximationVisualizer::calculateBatterPositionCorr(qint32 number, Mat img)
{

    Mat result;
    Mat pattern = cv::imread(QString("calibrate/%1batter_position.png").arg(number).toStdString());
    Mat searchZone;
    img(batterSearchZones[number]).copyTo(searchZone);
    qint32 result_cols =  searchZone.cols - pattern.cols + 1;
    qint32 result_rows = searchZone.rows - pattern.rows + 1;

    result.create( result_rows, result_cols, CV_32FC1 );

    matchTemplate( searchZone, pattern, result, CV_TM_CCORR_NORMED);
    double minVal; double maxVal; Point minLoc; Point maxLoc;

    minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );
    return maxVal;

}

void ApproximationVisualizer::setResultPicture(Mat m)
{
    QMutexLocker locker(&scenarioMutex);
    scenarioData.resultPicture = m;
}

void ApproximationVisualizer::handleNewScenarioState(RepeatVisualizeState state, cv::Mat mat, double time, bool clear = true)
{
    //resetRepeatVisualise();
    QMutexLocker locker(&scenarioMutex);
    if (clear)
    {
        scenarioData.repeats.clear();
    }

    if (state == ShowInit)
    {
        scenarioData.initPicture = mat;
    }
    else if (state == ShowResultPicture)
    {
        scenarioData.resultPicture = mat;
    }

    currentRepeatState = state;
    if (!qFuzzyCompare(time, -1))
    {
        scenarioData.waitTime = time;
        scenarioData.timer.start();
    }
}

void ApproximationVisualizer::setCurrentRepeatState(RepeatVisualizeState state, double time, cv::Mat mat)
{
    static RepeatVisualizeState tmpCurrentState = Invalid;
    if (state == tmpCurrentState)
    {
        return;
    }
    tmpCurrentState = state;
    qDebug() << "STATE CHANGE" << state << time << scenarioData.timer.elapsed();
    if (qFuzzyCompare(scenarioData.waitTime, -1) || scenarioData.waitTime < scenarioData.timer.elapsed())
    {
        handleNewScenarioState(state, mat, time);
    }
    else
    {
        qint32 delta = scenarioData.waitTime - scenarioData.timer.elapsed();
        if (delta > 0)
        {
            if (state == ShowRepeat)
            {
                QMutexLocker locker(&scenarioMutex);
                scenarioData.repeats.clear();
                locker.unlock();
            }
            qDebug() << "SET DELAY";
            QTimer::singleShot(delta, this, [this, time, mat, state]()
            {
                qDebug() << "HANDLE DELAY";
                handleNewScenarioState(state, mat, time, false);
            });
        }
    }

}

void ApproximationVisualizer::setLastApproximation(QSharedPointer <BallApproximator> lastApprox)
{
    scenarioData.approx = lastApprox;
}

bool ApproximationVisualizer::appendRepeatFrame(qint32 camNum, Mat frame, double time)
{
    QMutexLocker lock(&scenarioMutex); // вот тут так себе, потому что задерживает приём повтора, надо подумать
    if (scenarioData.repeats.contains(camNum))
    {
        scenarioData.repeats[camNum].second.append(qMakePair(frame, time));
    }
    else

    {
        qDebug() << "REMOVE" << camNum;
        return false;
        //scenarioData.repeats.insert(camNum, qMakePair(-1, QLinkedList <QPair <cv::Mat, double>>()));
    }
    return true;
    //scenarioData.repeat.append(qMakePair(frame, time));
}



void ApproximationVisualizer::setRepeatMainInfo(qint32 camNum, double coef, double sTime)
{
    QMutexLocker lock(&scenarioMutex);
    if (scenarioData.repeats.contains(camNum))
    {
        scenarioData.repeats.remove(camNum);
    }
    qDebug() << "SET CORR COEFF" << camNum << coef;
    scenarioData.repeats.insert(camNum, qMakePair(CorrTime(coef, -1, sTime), QLinkedList <QPair <cv::Mat, double>>()));
}

void ApproximationVisualizer::setRepeatFrameCount(qint32 camNum, qint32 count)
{
    if (scenarioData.repeats.contains(camNum))
    {
        scenarioData.repeats[camNum].first.frameCount = count;
    }
    else
    {
        Q_ASSERT(false);
    }
}

void ApproximationVisualizer::setRepeatCameraNumber(qint32 number)
{
    scenarioData.cameraNumber = number;
    CalibrationAdjustHelper::readCurrentCalibrationParameters(scenarioData.cameraNumber, "calibrate", scenarioData.EOcamera, scenarioData.camera, true);
}

void ApproximationVisualizer::setRepeatInitTime(qint32 num, double time)
{
    if (scenarioData.repeats.contains(num))
    {
        scenarioData.repeats[num].first.initTime = time;
    }
    else
    {
        Q_ASSERT(false);
    }
}

void ApproximationVisualizer::clearRepeats()
{
    QMutexLocker lock(&scenarioMutex);
    scenarioData.repeats.clear();
}

qint32 ApproximationVisualizer::correctCoordinates(double coord)
{
    return (coord - 956) * 1.33 + 1300;
}

void ApproximationVisualizer::correctPixmap(QPixmap &px, double &xOld)
{
    double xTmp = xOld;
    xOld = correctCoordinates(xOld);
    qint32 xScale = correctCoordinates(xTmp + px.width());
    px = px.scaled(xScale - xOld, px.height());
}


