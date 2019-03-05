#include "qchartviewwav.h"
#include "qwidgetwavchart.h"
#include <QImage>
#include "../../fftw-3.3.5-dll32/fftw3.h"

/*!
    \class QWidgetWAVChart
    The QWidgetWAVChart class inherits QWidget to show data.
*/
QWidgetWAVChart::QWidgetWAVChart(QWidget* parent) : QWidget(parent) {

    setFocusPolicy(Qt::StrongFocus);

    parentwindow = (MainWindow*)parentWidget();
    appbusy = false;
    nPoints = 600;
    detectCirclesize = 6;
    nadditionalfeatures = 0;
    hitii = -1;
    historyind = 0;
    Nseries = 0;

    initParams();

}

//QWidgetWAVChart::~QWidgetWAVChart(){};

void QWidgetWAVChart::initParams(){

    // *** Initialize figure parameters
    chartLayout = new QGridLayout;
    //chartLayout->setRowMinimumHeight(0,20);  // not working
    //chartLayout->setContentsMargins(0,0,0,0);  // not working
    //chartLayout->setSpacing(0);  // not working
    //chartLayout->setSizeConstraint(QLayout::SetFixedSize);

    seriespen.resize(1);
    seriespen[0].setWidth(1);
    seriespen[0].setColor(QColor(60,60,60,150));
    axisX = new QValueAxis();

    // *** Initialize figures
    axesYseries.resize(1);
    chart = new QChart();
    chart->legend()->hide();
    chart->setMargins(QMargins(0,0,0,0));
    //chart->setMinimumHeight(20.0); // with this, the charts sometimes alter the offset weirdly
    //QSizeF am = chart->sizeHint(Qt::MinimumSize);
    //chart->setMinimumSize(am);
    chart->setAnimationOptions(QChart::SeriesAnimations);
    //chart->setAnimationOptions(QChart::NoAnimation);

    //axesX.append(new QValueAxis);   //QObject::connect(chart, SIGNAL(scaleChanged()), this, SLOT(xaxisChanged()));
    chart->addAxis(axisX, Qt::AlignBottom);
    axesY.append(new QValueAxis);
    axisYspec.append(new QValueAxis);
    axisYspec[0]->setRange(0, 20000);
    chart->addAxis(axesY[0], Qt::AlignLeft);
    chart->addAxis(axisYspec[0], Qt::AlignRight);

    chartview = new QChartViewWAV(chart);
    chartview->setRenderHint(QPainter::Antialiasing);
    chartview->setID(0);

    //chart->setSizePolicy(QSizePolicy::ShrinkFlag);
    chart->setPreferredHeight(200); // with this, the charts sometimes alter the offset weirdly
    //chart->setPreferredWidth(1000); // with this, the charts sometimes alter the offset weirdly
    //chart[1]->setPreferredHeight(170); // with this, the charts sometimes alter the offset weirdly
    //chart[1]->setPreferredWidth(1000); // with this, the charts sometimes alter the offset weirdly

    // display default chart
        //QMessageBox::information(this, "!!", QString::number(nChart));
    chartLayout->addWidget(chartview, 0, 0);

    this->setLayout(chartLayout);
    axesYmin.resize(1);
    axesYmax.resize(1);

    axesY[0]->setVisible(false);
    axesYmin[0] = -30000;
    axesYmax[0] = 10000;
}


void QWidgetWAVChart::setData(){
    // MainWindow opened a wave file
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    QVector<double>& wavt = parentwindow->getDatat();
    Nseries = 1;

    duron.clear();
    duroff.clear();

    // delete all series and remove additional Yaxes
    for (int ii=series.length()-1; ii>=0; --ii){
        chart->removeSeries(series[ii]);
        delete(series.takeAt(ii));
    }
    if (axesY.length() > 1){
        for (int ii=1; ii<axesY.length(); ++ii){
            chart->removeAxis(axesY[ii]);
            axesY[ii]->deleteLater();
        }
        axesY.resize(1);
        axesYmin.resize(1);
        axesYmax.resize(1);
    }
    axesYseries.clear();
    axesYseries.resize(1);
    axesYseries[0].append(0); // series[0], which is the raw wav
    axesYmin[0] = -30000;
    axesYmax[0] = 10000;
    specoffsetx=0;
    specnumx=0;
    specstindx=0;
    specoffsety=0;
    specnumy=0;
    specstindy=0;

    // Add series[0]
    series.append(new QLineSeries());
    series[0]->setPen(seriespen[0]);   //series->setUseOpenGL(true);
    chart->addSeries(series[0]);
    //series[ii]->attachAxis(axesX[ii]);
    series[0]->attachAxis(axisX);
    series[0]->attachAxis(axesY[0]);

    axisX->setRange(wavt[0],wavt[wavt.length()-1]);

    updateSeriesX();
    saveAxesXHistory();

    //calcSpectrogram();
}


void QWidgetWAVChart::setSeries(int wavind, int axesYind){
    Nseries = wavind+1;
    //QMessageBox::information(this, "Nseries", QString::number(Nseries)); //2
    //QMessageBox::information(this, "axesY len", QString::number(axesY.length()-1)); //0

    if (axesYind > axesY.length()-1){
        axesY.append(new QValueAxis);
        axesYmin.append(0);
        axesYmax.append(0);
        chart->addAxis(axesY[axesYind], Qt::AlignLeft);
        axesYseries.resize(axesYind+1);
    }
    axesYseries[axesYind].append(Nseries-1);

    series.append(new QLineSeries());
    seriespen.resize(Nseries);
    seriespen[Nseries-1].setWidth(1);
    seriespen[Nseries-1].setColor(QColor(200,60,60,150));
    series[Nseries-1]->setPen(seriespen[Nseries-1]);   //series->setUseOpenGL(true);
    chart->addSeries(series[Nseries-1]);
    series[Nseries-1]->attachAxis(axisX);
    series[Nseries-1]->attachAxis(axesY[axesYind]);

    updateSeriesX();
}


void QWidgetWAVChart::saveData(QString fileName){
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    //QString rhdfilename = parentwindow->getFileName();
    QString rhdfilename = parentwindow->getFileName();
    QString rhdpathname = parentwindow->getFilePathName();
    rhdfilename.chop(4); // remove ".rhd"

    if (fileName.isEmpty()){
        fileName = QFileDialog::getSaveFileName(this, "Save as:", rhdpathname + "/" + rhdfilename + "_data.txt", "ASCII (*.txt);;All Files (*)");
    }
    rhdfilename= rhdfilename + ".rhd";
    if (fileName.isEmpty()){
        return;
    } else {
        QFile file(fileName);
        if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, "Unable to save file", file.errorString());
            return;
        }

        // write parametes
        QTextStream out(&file);
        out << "saveData" << "\r\n";
        out << rhdfilename << "\r\n";
        out << rhdpathname << "\r\n";
        // save something
        for(int ii=0; ii<Spec_phase.length(); ii++){
            for(int jj=0; jj<Spec_phase[ii].length(); jj++){
                out << Spec_phase[ii][jj]; // real value
                if (jj<Spec_phase[ii].length()-1){
                    out << "\t";
                }
            }
            out << "\r\n";
        }

        file.close();
    }
}



void QWidgetWAVChart::updateSeriesX(){

    // *** draw plots (series) in the figures
    if (!parentwindow->getFileOpened()){
        return;
    }

    int stpnt;
    int enpnt;
    int wavtlen;
    for (int ii=0; ii<Nseries; ++ii){
        series[ii]->clear();
    }

    QVector<double>& wavt = parentwindow->getDatat();

    qreal minaxisX = axisX->min();
    qreal maxaxisX = axisX->max();

    int minpnt = int(((minaxisX - wavt[0]) / (wavt[1] - wavt[0])));
    int maxpnt = int(((maxaxisX - wavt[0]) / (wavt[1] - wavt[0])));

    double minx;
    double miny;
    double maxx;
    double maxy;
    if (minpnt < 0){
        stpnt = 0;
    } else {
        stpnt = minpnt;
    }
    if (maxpnt + 1 > wavt.length()){
        enpnt = wavt.length();
    } else {
        enpnt = maxpnt + 1;
    }
    wavtlen = enpnt - stpnt;

    //QMessageBox::information(this, "!!tlength", "len:" + QString::number(wavtlen) + " range:" + QString::number(stpnt) + "-" + QString::number(enpnt));
    //QMessageBox::information(this, "!!tlength", "x:" + QString::number(wavt[stpnt]) + "-" + QString::number(datay[electrodearray[ii+stChart]][stpnt]) + " y:" + QString::number(miny) + "-" + QString::number(maxy));

    for (int ii=0; ii<axesY.length(); ++ii){
    // for (int ii=stChart; ii<nChart+stChart; ++ii){
    // *** ii = #Ch
        //QMessageBox::information(this, "!!ii:", QString::number(ii));
        //qreal minaxisX = axesX[ii]->min();
        //qreal maxaxisX = axesX[ii]->max();
        double allminy;
        double allmaxy;
        for (int kk=0; kk<axesYseries[ii].length(); ++kk){
            //QMessageBox::information(this, "!!kk:", QString::number(kk));
            //QMessageBox::information(this, "!!axesYseries[ii][kk]:", QString::number(axesYseries[ii][kk]));
            QVector<double>& anawav = parentwindow->getSeries(axesYseries[ii][kk]);
            miny = anawav[stpnt];
            maxy = anawav[stpnt];
            minx = wavt[stpnt];
            maxx = wavt[stpnt];
            allminy = miny;
            allmaxy = maxy;

            if (wavtlen <= nPoints){
                for (int jj = stpnt; jj < enpnt; ++jj){
                    series[axesYseries[ii][kk]]->append(wavt[jj], anawav[jj]);
                    //seriesfilt[axesYseries[ii][kk]]->append(wavt[jj], datayfilt[electrodearray[ii]][jj]);
                    if (allminy > anawav[jj]){
                        allminy = int(anawav[jj]);
                    }
                    if (allmaxy < anawav[jj]){
                        allmaxy = int(anawav[jj]);
                    }
                }
            } else {
                for (int jj = stpnt; jj < enpnt; ++jj){
                    if (jj%int(wavtlen/nPoints)==0){
                        if (minx < maxx){
                            series[axesYseries[ii][kk]]->append(minx, miny);
                            series[axesYseries[ii][kk]]->append(maxx, maxy);
                            //seriesfilt[axesYseries[ii][kk]]->append(minx, minyfilt);
                            //seriesfilt[axesYseries[ii][kk]]->append(maxx, maxyfilt);
                        } else {
                            series[axesYseries[ii][kk]]->append(maxx, maxy);
                            series[axesYseries[ii][kk]]->append(minx, miny);
                            //seriesfilt[axesYseries[ii][kk]]->append(maxx, maxyfilt);
                            //seriesfilt[axesYseries[ii][kk]]->append(minx, minyfilt);
                        }
                        if (allminy > miny){
                            allminy = (int)miny;
                        }
                        if (allmaxy < maxy){
                            allmaxy = (int)maxy;
                        }
                        miny = anawav[jj];
                        maxy = anawav[jj];
                        minx = wavt[jj];
                        maxx = wavt[jj];
                        //maxyfilt = datayfilt[electrodearray[axesYseries[ii][kk]]][jj];
                        //minyfilt = datayfilt[electrodearray[axesYseries[ii][kk]]][jj];
                    } else {
                        if (maxy < anawav[jj]){
                            maxy = anawav[jj];
                            maxx = wavt[jj];
                            //maxyfilt = datayfilt[electrodearray[axesYseries[ii][kk]]][jj];
                        }
                        if (miny > anawav[jj]){
                            miny = anawav[jj];
                            minx = wavt[jj];
                            //minyfilt = datayfilt[electrodearray[axesYseries[ii][kk]]][jj];
                        }
                    }
                }
            }
        }

        if (axesYmin[ii] == 0 && axesYmax[ii] == 0){
            axesY[ii]->setRange((allminy - int(allminy)%100 - 100), (allmaxy - int(allmaxy)%100 + 100));
        } else {
            axesY[ii]->setRange(axesYmin[ii], axesYmax[ii]);
        }


        //drawDetect(ii+stChart);


        //axesX[ii]->setRange(mintick, maxtick);
        //axesX[ii]->setTickCount(ntick);

        //if (seriesthres[ii]->count() > 0){
        //    seriesthres[indCh]->clear();
        //    seriesthres[indCh]->append(axesX[indCh]->min(), thres);
        //    seriesthres[indCh]->append(axesX[indCh]->max(), thres);
        //    seriesdetect[indCh]->clear();
        //}
    }

    /* test
    if (series.length()>1){
        QVector<double> vec;
        vec.resize(series[1]->count());
        for (int ii=0; ii<series[1]->count(); ++ii){
            vec[ii] = series[1]->at(ii).x();
        }
        parentwindow->saveData1D(vec);
    }
    */

    //QMessageBox::information(this, "!!tlength", QString::number(wavt.length()));
    //QMessageBox::information(this, "!!stpnt", QString::number(stpnt));
    //QMessageBox::information(this, "!!enpnt", QString::number(enpnt));

    prepareSpecImage();
}

void QWidgetWAVChart::setAxesX(qreal min, qreal max){
    QVector<double>& wavt = parentwindow->getDatat();
    if (min == 0 && max == 0){
        min = wavt[0];
        max = wavt[wavt.length() - 1];
    }
    if (min >= max){
        return;
    }
    //for (int ii=0; ii<Nseries; ++ii){
        //axesX[ii]->setRange(min,max);
    //}
    axisX->setRange(min,max);
    updateSeriesX();
}

void QWidgetWAVChart::setAxisYspec(int ind, qreal min, qreal max){
    axisYspec[ind]->setRange(min, max);
}

void QWidgetWAVChart::setAxesY(){
    QVector<QString> dialoglabel;
    QVector<QString> dialogvalue;
    QVector<QString> returnvalue;
    dialoglabel.append("axis num (0 for all)");
    dialoglabel.append("Y max");
    dialoglabel.append("Y min (reset if max&min == 0)");
    dialogvalue.append("0");
    /*
    dialogvalue.append(QString::number(axesYmax[0]));
    dialogvalue.append(QString::number(axesYmin[0]));
    returnvalue = myDialog("Y range:", dialoglabel, dialogvalue);
    if (returnvalue.length() == 0){
        return;
    }
    int axisnum = returnvalue[0].toInt();
    qreal ymax = returnvalue[1].toDouble();
    qreal ymin = returnvalue[2].toDouble();
    if (axisnum == 0){
        for (int ii=0; ii<nChart; ++ii){
            axesYmax[ii] = ymax;
            axesYmin[ii] = ymin;
        }
    } else if (axisnum >= 1 && axisnum <= 16) {
        axesYmax[axisnum-1] = ymax;
        axesYmin[axisnum-1] = ymin;
    }
    */
    updateSeriesX();
}

void QWidgetWAVChart::scaleAxesX(double scalemin, double scalemax){
    /*
    for (int ii=0; ii<Nseries; ++ii){
        qreal min = axesX[ii]->min();
        qreal max = axesX[ii]->max();
        qreal len = max - min;
        axesX[ii]->setRange(min + len*scalemin, max + len*scalemax);
    }
    */
    qreal min = axisX->min();
    qreal max = axisX->max();
    qreal len = max - min;
    axisX->setRange(min + len*scalemin, max + len*scalemax);
    updateSeriesX();
}

void QWidgetWAVChart::alignAxesXto(int indChart){
    /*
    qreal min = axesX[indChart]->min();
    qreal max = axesX[indChart]->max();
    for (int ii=0; ii<Nseries; ++ii){
        axesX[ii]->setRange(min,max);
    }
    updateSeriesX();
    saveAxesXHistory();
    */
    updateSeriesX();
    saveAxesXHistory();
}

void QWidgetWAVChart::saveAxesXHistory(){
    //QMessageBox::information(this, "!prehistoryind", QString::number(historyind));
    if (historyind == 10){
        for (int ii=0; ii<9; ++ii){
            historymin[ii] = historymin[ii+1];
            historymax[ii] = historymax[ii+1];
        }
        historymin[9] = axisX->min();
        historymax[9] = axisX->max();
    } else {
        if (historyind == historymin.size()){
            historymin.append(axisX->min());
            historymax.append(axisX->max());
            historyind += 1;
        } else {
            historymin[historyind] = axisX->min();
            historymax[historyind] = axisX->max();
            historyind += 1;
        }
    }
    //QMessageBox::information(this, "!posthistoryind", QString::number(historyind));
}


void QWidgetWAVChart::keyPressReceiver(QKeyEvent* event){
    switch (event->key()) {
        case Qt::Key_Home:
            setAxesX(0, 0);
            saveAxesXHistory();
            break;
        case Qt::Key_0:
            if (event->modifiers() == Qt::ShiftModifier){
                setAxesX(-0.01, 0.03);
            } else {
                setAxesX(-0.02, 0.06);
            }
            saveAxesXHistory();
            break;
        case Qt::Key_Up:
            if (event->modifiers() == Qt::ShiftModifier){
                scaleAxesX(3.5/8.0, -3.5/8.0);
            } else {
                scaleAxesX(1.0/4.0, -1.0/4.0);
            }
            saveAxesXHistory();
            break;
        case Qt::Key_Down:
            if (event->modifiers() == Qt::ShiftModifier){
                scaleAxesX(-3.5, 3.5);
            } else {
                scaleAxesX(-1.0/2.0, 1.0/2.0);
            }
            saveAxesXHistory();
            break;
        case Qt::Key_Left:
            if (event->modifiers() == Qt::ShiftModifier){
                scaleAxesX(-1.0/2.0, -1.0/2.0);
            } else {
                scaleAxesX(-1.0/8.0, -1.0/8.0);
            }
            saveAxesXHistory();
            break;
        case Qt::Key_Right:
            if (event->modifiers() == Qt::ShiftModifier){
                scaleAxesX(1.0/2.0, 1.0/2.0);
            } else {
                scaleAxesX(1.0/8.0, 1.0/8.0);
            }
            saveAxesXHistory();
            break;
        case Qt::Key_Z:
            if (event->modifiers() == Qt::ControlModifier){
                if (historyind > 1){
                    //QMessageBox::information(this, "!prehistoryind", QString::number(historyind));
                    historyind -= 1;
                    setAxesX(historymin[historyind-1], historymax[historyind-1]);
                    //QMessageBox::information(this, "!posthistoryind", QString::number(historyind));
                }
            }
            break;
        case Qt::Key_Y:
            if (event->modifiers() == Qt::ControlModifier){
                if (historyind < historymin.size()){
                    //QMessageBox::information(this, "!prehistoryind", QString::number(historyind));
                    historyind += 1;
                    setAxesX(historymin[historyind-1], historymax[historyind-1]);
                    //QMessageBox::information(this, "!posthistoryind", QString::number(historyind));
                }
            }
            break;
        case Qt::Key_D:
            if (event->modifiers() == Qt::ControlModifier){
                setAxesY();
            }
            break;
    }
}

void QWidgetWAVChart::prepareSpecImage(){
    // still the spectrogram is a bit offseted (when the number of windows is large. I think this is due to the rounding errors (probably dx is not accurate))

    //struct spec& sp = parentwindow->getSpec();
    Spec_mag_max = parentwindow->getSpec_mag_max();
    QVector<double>& Spec_x = parentwindow->getSpec_x();
    QVector<double>& Spec_y = parentwindow->getSpec_y();
    Spec_dx = parentwindow->getSpec_dx();
    Spec_dy = parentwindow->getSpec_dy();

    if (Spec_x.length() == 0) return;
    //QMessageBox::information(this, "!prepareSpecImage", "okok");
    //axisYspec[0]->setRange(0, 20000);
    // left & right t
    qreal minx = axisX->min(); // mint ~ Topleft.x()
    qreal maxx = axisX->max(); // maxt ~ Topleft.x() + width
    qreal miny = axisYspec[0]->min(); // mint ~ Topleft.x()
    qreal maxy = axisYspec[0]->max(); // maxt ~ Topleft.x() + width
    double xwid = maxx - minx;
    double ywid = maxy - miny;

    specoffsetx=0;
    specnumx=0;
    specstindx=0;
    for (int ii=0; ii<Spec_x.length(); ii++){
        if (minx <= Spec_x[ii] && Spec_x[ii] <= maxx){
            if (specnumx == 0){
                specoffsetx = Spec_x[ii] - minx;
                specstindx = ii;
            }
            specnumx++;
        } else if (Spec_x[ii] > maxx){
            break;
        }
    }
    //QMessageBox::information(this, "!specnumx", QString::number(specnumx));
    if (specnumx == 0) return;
    specoffsety=0;
    specnumy=0;
    specstindy=0;

    for (int ii=0; ii<Spec_y.length(); ii++){
        if (miny <= Spec_y[ii] && Spec_y[ii] <= maxy){
            if (specnumy == 0){
                specoffsety = Spec_y[ii] - miny;
                specstindy = ii;
            }
            specnumy++;
        } else if (Spec_y[ii] > maxy){
            break;
        }
    }
    //QMessageBox::information(this, "!posthistoryind", QString::number(specnumy));
    if (specnumy == 0) return;

    specoffsetx /= xwid;
    specoffsety /= ywid;
    specboxx = Spec_dx/xwid;
    specboxy = Spec_dy/ywid;

    //QMessageBox::information(this, "!specoffsetx", QString::number(specoffsetx));
    //QMessageBox::information(this, "!specoffsety", QString::number(specoffsety));
    //QMessageBox::information(this, "!specboxx", QString::number(specboxx));
    //QMessageBox::information(this, "!specboxy", QString::number(specboxy));

    // *** prepareDuration
    QVector<double>& ton = parentwindow->getONt();
    QVector<double>& toff = parentwindow->getOFFt();
    //QMessageBox::information(this, "!specboxy", QString::number(ton.length()));

    if (ton.length() != 0 && toff.length() != 0){
        duron.clear();
        duroff.clear();
        for (int ii=0; ii<ton.length(); ii++){
            if ((minx <= ton[ii] && ton[ii] <= maxx) || (minx <= toff[ii] && toff[ii] <= maxx)){
                duron.append((ton[ii]-minx)/xwid);
                duroff.append((toff[ii]-minx)/xwid);
                //QMessageBox::information(this, "!data", "[" + QString::number(ii) + "]\r\n" + QString::number(duron[ii]) + "-" + QString::number(duroff[ii]));
            } else if (ton[ii] > maxx){
                break;
            }
        }
    }

    repaint();
}



void QWidgetWAVChart::paintEvent(QPaintEvent* event){
    // Remember that paintEvent() can be called anytime. This often generates errors.
    QVector<QVector<double>>& Spec_mag = parentwindow->getSpec_mag();
    if (Spec_mag.length() == 0 || specnumx == 0 || specnumy == 0) { // in case Spec_mag is cleared in mainwindow.cpp
        return;
    } else if (parentwindow->appbusy){
        return;
    }

    // prepare image space for spectrogram
    qreal plotwidth = chart->plotArea().width();
    qreal plotheight = chart->plotArea().height();
    qreal viewwidth = chartview->width();
    qreal viewheight = chartview->height();    // Spectrograms

    QImage imgspec(viewwidth, viewheight, QImage::Format_ARGB32); // rather than plotwidth
    QPointF TopLeft = chart->plotArea().topLeft();
    QPointF BottomLeft = chart->plotArea().bottomLeft();
    //painter.drawImage(TopLeft, imgspec);
    //painter.drawImage(50,50, imgspec);

    // paint spectrogram (this is not optimal because the calc is more than resolution)
    imgspec.fill(Qt::white);
    QPainter mypainter(&imgspec);
    QRect tmprec;
    int gainshade = parentwindow->getSpecslider();
    int colorval = 0;

    for (int ii=0; ii<specnumx; ii++){
        tmprec.setLeft(int(TopLeft.x() + (specoffsetx - specboxx/2 + specboxx * ii)*plotwidth));
        //tmprec.setLeft(int((specoffsetx - specboxx/2 + specboxx * ii)*plotwidth));
        tmprec.setWidth(int(specboxx * plotwidth)+1);
        if (tmprec.width() == 0){
            tmprec.setWidth(1);
        }
        for (int jj=0; jj<specnumy; jj++){
            tmprec.setTop(int(TopLeft.y() + (1 - (specoffsety + specboxy/2 + specboxy * jj))*plotheight));
            //tmprec.setTop(int((1 - (specoffsety + specboxy/2 + specboxy * jj))*plotheight));
            tmprec.setHeight(int(specboxy * plotheight)+1);
            if (tmprec.height() == 0){
                tmprec.setHeight(1);
            }
            colorval = 255 - int(Spec_mag[ii+specstindx][jj+specstindy]/Spec_mag_max * 255) * gainshade;
            if (colorval < 0) colorval = 0;
            mypainter.fillRect(tmprec, QColor(colorval, colorval, colorval));
//            mypainter.fillRect(tmprec, Qt::red);
        }
    }

    // drawDuration
    for (int ii=0; ii<duron.length(); ii++){
        tmprec.setLeft(int(BottomLeft.x() + (duron[ii])*double(plotwidth)));
        tmprec.setTop(int(BottomLeft.y() - 20));
        tmprec.setHeight(20);
        tmprec.setWidth(int((duroff[ii] - duron[ii])*double(plotwidth)));
        mypainter.fillRect(tmprec, QColor(0, 0, 200, 50));
    }

    /*
    tmprec.setLeft(int(TopLeft.x()));
    //tmprec.setLeft(0);
    tmprec.setTop(int(TopLeft.y()));
    //tmprec.setTop(0);
    tmprec.setWidth(int(plotwidth));
    tmprec.setHeight(int(plotheight));
    mypainter.fillRect(tmprec, Qt::red);
    */

    QBrush imgbrush(imgspec);
    chart->setPlotAreaBackgroundBrush(imgbrush);
    chart->setPlotAreaBackgroundVisible(true);


    /*
    QRgb value = qRgb(0, 200, 0);
    imgspec.setColor(0, value);
    for (int ii=0; ii<50; ii++){
        for (int jj=0; jj<50; jj++){
            imgspec.setPixel(ii,jj,0);
        }
    }
    */

    //chart[1]->legend()->hide();
    //chart[1]->createDefaultAxes();
    //chart[1]->axisX()->setGridLineVisible(false);
    //chart[1]->axisY()->setGridLineVisible(false);
}

QVector<QVector<int>> QWidgetWAVChart::getAxesYseries(){
    return axesYseries;
}

qreal QWidgetWAVChart::getXmax(){
    return axisX->max();
}

qreal QWidgetWAVChart::getXmin(){
    return axisX->min();
}

void QWidgetWAVChart::reprint(){
    repaint();
}
