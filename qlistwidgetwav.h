#ifndef QLISTWIDGETWAV_H
#define QLISTWIDGETWAV_H

#include "mainwindow.h"

class QListWidgetWAV : public QListWidget {
    public:
        QListWidgetWAV(QWidget* = nullptr);
    private:
        MainWindow* parentwindow;
};

#endif // QLISTWIDGETWAV_H
