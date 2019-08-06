#ifndef EDITFRAMEQGRAPHICSVIEW_H
#define EDITFRAMEQGRAPHICSVIEW_H

#include <QObject>
#include <QWidget>
#include <QGraphicsView>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsRectItem>
#include <QWheelEvent>
#include <QtMath>

class EditFrameQGraphicsScene : public QGraphicsScene
{
    Q_OBJECT
public:
    EditFrameQGraphicsScene(QObject* parent = nullptr);
    ~EditFrameQGraphicsScene();
    void setFrame(const QPixmap& frame)
    {
        mainFrame->setPixmap(frame);
    }

    void clearROI();

    QRect getROI() {return rect;}

    void drawROI(QRect dRect) {drawROIInternal(dRect);}
signals:
    void ROISelected(QRect rect);

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);

    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);

    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

    void keyPressEvent(QKeyEvent* event);

    void wheelEvent(QWheelEvent* pWheelEvent);

private:
    QPointF rectStart;
    QPointF rectEnd;
    QRect rect;
    QGraphicsPixmapItem* mainFrame = nullptr;
    void drawROIInternal(QRect dRect);
};

#endif // EDITFRAMEQGRAPHICSVIEW_H
