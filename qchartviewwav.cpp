#include "qchartviewwav.h"

/*!
    \class QChartViewWAV
    \brief The QChartViewWAV class inherits QChartView to show data.
    \since 0.0
    \ingroup main class

    \sa MainWindow, QWidgetWAVChart
*/

QChartViewWAV::QChartViewWAV (QChart *chart, QWidget *parent) : QChartView(chart, parent), m_isTouching(false){
    setRubberBand(QChartView::HorizontalRubberBand);
    mouserelease_busy = false;
    ID = -1;
    parentwidget = (QWidgetWAVChart*) parentWidget();
}

void QChartViewWAV::setID(int tmpID){
    ID = tmpID;
}

int QChartViewWAV::getID(){
    return ID;
}

void QChartViewWAV::keyPressEvent(QKeyEvent *event){
    parentwidget = (QWidgetWAVChart*) parentWidget();
    parentwidget->keyPressReceiver(event);
}


bool QChartViewWAV::viewportEvent(QEvent *event){
    if (event->type() == QEvent::TouchBegin) {
        // By default touch events are converted to mouse events. So
        // after this event we will get a mouse event also but we want
        // to handle touch events as gestures only. So we need this safeguard
        // to block mouse events that are actually generated from touch.
        m_isTouching = true;

        // Turn off animations when handling gestures they
        // will only slow us down.
        chart()->setAnimationOptions(QChart::NoAnimation);
    }
    return QChartView::viewportEvent(event);
}

void QChartViewWAV::mousePressEvent(QMouseEvent *event){
    if (m_isTouching){
        return;
    }
    pospressed = QWidget::mapFromGlobal(QCursor::pos());
    QChartView::mousePressEvent(event);
}

void QChartViewWAV::mouseMoveEvent(QMouseEvent *event){
    if (m_isTouching){
        return;
    }
    QChartView::mouseMoveEvent(event);
}

void QChartViewWAV::mouseReleaseEvent(QMouseEvent *event){
    if (m_isTouching){
        m_isTouching = false;
    }

    if (mouserelease_busy){
        return;
    } else {
        if (event->button() == Qt::RightButton){

        } else {
            QChartView::mouseReleaseEvent(event);
            mouserelease_busy = true; // This is to avoid computations when users click rapidly. I am not sure this works...

            //parentwidget = (WidgetIntanChart*)parentWidget();

            // Because we disabled animations when touch event was detected
            // we must put them back on.
            //chart()->setAnimationOptions(QChart::SeriesAnimations);

            QPointF pos = QWidget::mapFromGlobal(QCursor::pos());
            if (pos.x() == pospressed.x() && pos.y() == pospressed.y()){
                //QMessageBox::information(this, "hit!!", QString::number(pos.x()) + "," + QString::number(pos.y()));
                pos = chart()->mapFromParent(pos);
                //QMessageBox::information(this, "hit!!", QString::number(pos.x()) + "," + QString::number(pos.y()));
                pos = chart()->mapToValue(pos);
                //QMessageBox::information(this, "hit!!", QString::number(pos.x()) + "," + QString::number(pos.y()));


                //parentwidget->calcHitSpike(pos);

            } else {
                //flagmousepressed = false;
                parentwidget = (QWidgetWAVChart*) parentWidget(); // crash if this is omitted. QWidgetWAVChart* changes after entering the constructor
                parentwidget->alignAxesXto(ID);
            }
            mouserelease_busy = false;
        }
    }
}
