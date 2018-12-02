#ifndef QWIDGETWAVCHART_H
#define QWIDGETWAVCHART_H
#include <QWidget>
#include "mainwindow.h"
//#include "./qchartviewwav.h"

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
    private:
        MainWindow* parentwindow;
        QVector<QChartViewWAV*> chartview;
        QChartViewWAV* chartview0;

        // variables
        const static int nChart = 3;
        const static int stChart = 0;
        const static int maxChart = 3;
        const double PI = 3.141592653589793238463;
        bool appbusy;
        int nPoints;
        int detectCirclesize;
        int nadditionalfeatures;
        int hitii;
        int nChannel;
        QVector<qreal> historymin;
        QVector<qreal> historymax;
        int historyind = 0;

        // bunch of lists // QVector would be better
        QVector<QChart*> chart;
        QGridLayout *chartLayout;
        QVector<QLineSeries*> series;
        QVector<QValueAxis*> axesX;
        QVector<QValueAxis*> axesY;
        QPen linepen;
        QVector<QWidget*> buttons;
        QVector<QSignalMapper*> buttonsignalmapper; // used to be [16][3] array, but changed to a vector
        QVector<QPushButton*> buttonShowFilt;
        QVector<QPushButton*> buttonDetectSpike;
        QVector<QPushButton*> buttonShowSpike;
        QVector<QVBoxLayout*> buttonLayout;
        QVector<bool> buttonfiltpressed;

        QVector<qreal> axesYmin;
        QVector<qreal> axesYmax;

        QVector<QVector<double>> Spec_mag;
        double Spec_mag_max;
        QVector<QVector<double>> Spec_phase;
        QVector<double> Spec_x;
        QVector<double> Spec_y;
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
        void updateSeriesX();
        void setAxesX(qreal, qreal);
        void setAxesY();
        void scaleAxesX(double, double);
        void saveAxesXHistory();
        void initParams();
        void calcSpectrogram();
        void prepareSpecImage();
    protected:
        void paintEvent(QPaintEvent*);

};



#endif // QWIDGETWAVCHART_H
