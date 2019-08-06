#ifndef ZOOMGRAPHICSVIEW_H
#define ZOOMGRAPHICSVIEW_H

#include <QGraphicsView>
#include <QWheelEvent>

class ZoomGraphicsView : public QGraphicsView
{
    Q_OBJECT
public:
    ZoomGraphicsView(QWidget* parent = nullptr);
protected:
    void wheelEvent(QWheelEvent* pWheelEvent);
};

#endif // ZOOMGRAPHICSVIEW_H
