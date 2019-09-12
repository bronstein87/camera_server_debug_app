#include "editframeqgraphicsview.h"

EditFrameQGraphicsScene::EditFrameQGraphicsScene(QObject* parent) : QGraphicsScene(parent),
   mainFrame(new QGraphicsPixmapItem())
{
    addItem(mainFrame);
}

EditFrameQGraphicsScene::~EditFrameQGraphicsScene()
{
    delete mainFrame;
}

void EditFrameQGraphicsScene::clearROI()
{
    for (auto& i : mainFrame->childItems())
    {
        removeItem(i);
        delete i;
    }
    rect = QRect();
}



void EditFrameQGraphicsScene::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    rectStart = event->scenePos();
    QGraphicsScene::mousePressEvent(event);
}

void EditFrameQGraphicsScene::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    QGraphicsScene::mouseMoveEvent(event);
}


#include <QDebug>
void EditFrameQGraphicsScene::drawROIInternal(QRect dRect)
{
    auto itemToDraw = new QGraphicsRectItem(mainFrame);
    itemToDraw->setPen(QPen(Qt::red, 3, Qt::SolidLine));
    itemToDraw->setRect(dRect);
}

void EditFrameQGraphicsScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    rectEnd = event->scenePos();
    rect.setX(rectStart.x());
    rect.setY(rectStart.y());
    rect.setWidth(rectEnd.x() - rectStart.x());
    rect.setHeight(rectEnd.y() - rectStart.y());
    drawROIInternal(rect);
    emit ROISelected(rect);
    QGraphicsScene::mouseReleaseEvent(event);
}

void EditFrameQGraphicsScene::keyPressEvent(QKeyEvent* event)
{
    if(event->key() == Qt::Key_Delete)
            clearROI();
    else
        QGraphicsScene::keyPressEvent(event);
}




