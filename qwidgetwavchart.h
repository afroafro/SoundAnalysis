#ifndef QWIDGETWAVCHART_H
#define QWIDGETWAVCHART_H
#include <QWidget>
#include "mainwindow.h"

class QChartViewWAV;
using namespace QtCharts; // necessary for some compilers
//class MainWindow;

class QWidgetWAVChart : public QWidget {
    Q_OBJECT
    public:
        QWidgetWAVChart(QWidget *parent = nullptr);
        //virtual ~QWidgetWAVChart() = 0;
        void keyPressReceiver(QKeyEvent*);
        void setData();
        void alignAxesXto(int);
        void saveData(QString=nullptr);
        void setSeries(int, int);
        QVector<QVector<int>> getAxesYseries();
        qreal getXmax();
        qreal getXmin();
        void prepareSpecImage();
        void setAxisYspec(int ind, qreal min, qreal max);
        void reprint();
    private:
        MainWindow* parentwindow;

        // variables
        const double PI = 3.141592653589793238463;
        int Nseries;
        bool appbusy;
        int nPoints;
        int detectCirclesize;
        int nadditionalfeatures;
        int hitii;
        QVector<qreal> duron;
        QVector<qreal> duroff;
        QVector<qreal> historymin;
        QVector<qreal> historymax;
        int historyind = 0;

        // bunch of lists // QVector would be better
        QChart* chart;
        QChartViewWAV* chartview;
        QGridLayout *chartLayout;
        QVector<QLineSeries*> series;
        QValueAxis* axisX;
        QVector<QValueAxis*> axisYspec;
        QVector<QValueAxis*> axesY;
        QVector<QVector<int>> axesYseries;
        QVector<QPen> seriespen;
        QVector<bool> buttonfiltpressed;

        QVector<qreal> axesYmin;
        QVector<qreal> axesYmax;

        QVector<QVector<double>> Spec_mag;
        double Spec_mag_max;
        QVector<QVector<double>> Spec_phase;
        double Spec_dx;
        double Spec_dy;
        QVector<QVector<double>> FFTresult;
        QVector<double> FFTmag;
        QVector<double> FFTphase;

        double specoffsetx;
        double specoffsety;
        double specboxx;
        double specboxy;
        int specnumx=0;
        int specnumy=0;
        int specstindx;
        int specstindy;

        // functions
        void setAxesX(qreal, qreal);
        void setAxesY();
        void scaleAxesX(double, double);
        void saveAxesXHistory();
        void initParams();
        void calcSpectrogram();
        void updateSeriesX();
    protected:
        void paintEvent(QPaintEvent*);

};



#endif // QWIDGETWAVCHART_H
