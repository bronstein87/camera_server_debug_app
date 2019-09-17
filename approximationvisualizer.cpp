#include "approximationvisualizer.h"

ApproximationVisualizer::ApproximationVisualizer(QObject *parent) : QObject(parent), scenePicture (new QGraphicsScene()),
    fullPixmap("://resources/picture.png"), shortPixmap("://resources/picture3.png")
{

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

void ApproximationVisualizer::plotCamera(BallApproximator& approx, qint32 firstNumber, qint32 secondNumber, qint32 recCountF, qint32 recCountS, QDateTime dt)
{
    QVector <double> timeFirst = approx.getFirstTime();
    QVector <double> timeSecond = approx.getSecondTime();
    QVector <double> timeFirstInit = approx.getFirstTimeInit();
    QVector <double> timeSecondInit = approx.getSecondTimeInit();
    QVector <double> errorsFirst = approx.getFirstErrors();
    QVector <double> errorsSecond = approx.getSecondErrors();
    QVector <double> timeFirstDelta;
    QVector <double> timeSecondDelta;

    fillSkippedTime(errorsFirst, timeFirstInit, timeFirst);
    fillSkippedTime(errorsSecond, timeSecondInit, timeSecond);

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
    m.objCount = recCountF;
    measuresMap[firstNumber].append(m);

    m.error = errorsSecond;
    m.timeDelta = timeSecondDelta;
    m.objCount = recCountS;
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
                                      xBegin, rot, W);
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

    return image;
}

void ApproximationVisualizer::drawZones(Calibration::ExteriorOr& eOr, Calibration::SpacecraftPlatform::CAMERA::CameraParams& cam)
{
    double widthInit = 0;
    const double closeZoneY = 0;
    const double farZoneY = 0.43;
    const double width = 0.215;
    const double minHeight = 0.473;
    const double maxHeight = 1.045;//1.075;

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
    double xFoneRes = 156;
    double yFoneRes = 188;
    QGraphicsPixmapItem* resFoneItem =  new QGraphicsPixmapItem(picture);
    resFoneItem->setPixmap(resFone);
    resFoneItem->setPos(xFoneRes, yFoneRes);

    QPixmap pitch("://resources/pitch.png");
    double xTitle = 250;
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
    auto effect =new QGraphicsDropShadowEffect(textItemVBegin);
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
                                                    xBegin, rot, W);
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
    return image;
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

