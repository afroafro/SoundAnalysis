#ifndef QCHARTVIEWWAV_H
#define QCHARTVIEWWAV_H

#include "qwidgetwavchart.h" // I have to include this because I want to use the member function of this class in .cpp
//#include <QtCharts>
//class QWidgetWAVChart;

class QChartViewWAV : public QChartView {
    public:
        QChartViewWAV(QChart *chart=nullptr, QWidget *parent=nullptr);
        void setID(int);
        int getID();
    private:
        bool m_isTouching;
        bool mouserelease_busy;
        int ID;
        QPointF pospressed;
        QWidgetWAVChart* parentwidget;
        void keyPressEvent(QKeyEvent *event);
    protected:
        void mousePressEvent(QMouseEvent *event);
        void mouseMoveEvent(QMouseEvent *event);
        void mouseReleaseEvent(QMouseEvent *event);
        bool viewportEvent(QEvent *event);
};

#endif // QCHARTVIEWWAV_H
