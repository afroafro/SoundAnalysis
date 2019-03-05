#include "mainwindow.h"

/*!
    \class main
    \brief The main class executes MainWindow.
    \since 0.0
    \ingroup main class

    \sa MainWindow, QChartViewWAV, QWidgetWAVChart
*/

int main(int argc, char *argv[]) {
    QApplication a(argc, argv);
    MainWindow w;
    w.grabGesture(Qt::PanGesture);
    w.grabGesture(Qt::PinchGesture);
    w.resize(400, 300);
    w.show();
    return a.exec();
}
