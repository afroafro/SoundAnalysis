#include "qchartviewwav.h"
#include "qwidgetwavchart.h"
#include <QImage>
#include "../../fftw-3.3.5-dll32/fftw3.h"

QWidgetWAVChart::QWidgetWAVChart(QWidget* parent) : QWidget(parent) {

    setFocusPolicy(Qt::StrongFocus);

    parentwindow = (MainWindow*)parentWidget();
    appbusy = false;
    nPoints = 600;
    detectCirclesize = 6;
    nadditionalfeatures = 0;
    hitii = -1;
    historyind = 0;
    nChannel = 1;

    initParams();

}

//QWidgetWAVChart::~QWidgetWAVChart(){};

void QWidgetWAVChart::initParams(){

    // *** Initialize figure parameters
    chartLayout = new QGridLayout;
    chartLayout->setRowMinimumHeight(50,20);
    chartLayout->setContentsMargins(0,0,0,0);
    chartLayout->setSpacing(0);

    linepen.setWidth(1);
    linepen.setColor(QColor(60,60,60,150));

    // *** Initialize figures
    for (int ii=0; ii<nChart; ++ii){
        chart.append(new QChart());
        chart[ii]->legend()->hide();
        chart[ii]->setMargins(QMargins(0,0,0,0));
        chart[ii]->setMinimumHeight(20.0);
        chart[ii]->setAnimationOptions(QChart::SeriesAnimations);
        //chart[ii]->setAnimationOptions(QChart::NoAnimation);

        axesX.append(new QValueAxis);   //QObject::connect(chart[ii], SIGNAL(scaleChanged()), this, SLOT(xaxisChanged()));
        chart[ii]->addAxis(axesX[ii], Qt::AlignBottom);
        axesY.append(new QValueAxis);
        chart[ii]->addAxis(axesY[ii], Qt::AlignLeft);

        // Buttons
        buttons.append(new QWidget());
        buttonDetectSpike.append(new QPushButton("Detect" + QString::number(ii+1)));
        buttonShowFilt.append(new QPushButton("Show filter" + QString::number(ii+1)));
        buttonShowSpike.append(new QPushButton("Show spikes" + QString::number(ii+1)));
        buttonfiltpressed.append(false);

        buttons[ii]->setEnabled(false);
        buttonDetectSpike[ii]->setEnabled(false);
        buttonShowFilt[ii]->setEnabled(false);
        buttonShowSpike[ii]->setEnabled(false);

        // to define the button functions
        buttonsignalmapper.append(new QSignalMapper(buttons[ii]));
        QObject::connect(buttonDetectSpike[ii], SIGNAL(released()), buttonsignalmapper[ii*3], SLOT(map()));
        buttonsignalmapper.append(new QSignalMapper(buttons[ii]));
        QObject::connect(buttonShowFilt[ii], SIGNAL(released()), buttonsignalmapper[ii*3+1], SLOT(map()));
        buttonsignalmapper.append(new QSignalMapper(buttons[ii]));
        QObject::connect(buttonShowSpike[ii], SIGNAL(released()), buttonsignalmapper[ii*3+2], SLOT(map()));
        buttonsignalmapper[ii*3]->setMapping(buttonDetectSpike[ii], ii); // if clicked throw ii (Channel name)
        buttonsignalmapper[ii*3+1]->setMapping(buttonShowFilt[ii], ii);
        buttonsignalmapper[ii*3+2]->setMapping(buttonShowSpike[ii], ii);
        QObject::connect(buttonsignalmapper[ii*3], SIGNAL(mapped(int)), this, SLOT(buttonpushDetectSpike(int)));
        QObject::connect(buttonsignalmapper[ii*3+1], SIGNAL(mapped(int)), this, SLOT(buttonpushShowFilt(int)));
        QObject::connect(buttonsignalmapper[ii*3+2], SIGNAL(mapped(int)), this, SLOT(buttonpushShowSpike(int)));

        buttonLayout.append(new QVBoxLayout);
        buttonLayout[ii] -> addWidget(buttonDetectSpike[ii]);
        buttonLayout[ii] -> addWidget(buttonShowFilt[ii]);
        buttonLayout[ii] -> addWidget(buttonShowSpike[ii]);
        buttons[ii] -> setLayout(buttonLayout[ii]);

        series.append(new QLineSeries());
        series[ii]->setPen(linepen);   //series->setUseOpenGL(true);
        chart[ii]->addSeries(series[ii]);
        series[ii]->attachAxis(axesX[ii]);
        series[ii]->attachAxis(axesY[ii]);

        chartview.append(new QChartViewWAV(chart[ii]));
        chartview[ii]->setRenderHint(QPainter::Antialiasing);
        chartview[ii]->setID(ii);
    }

    // display default chart
    for (int ii=0; ii<nChart; ++ii){
        //QMessageBox::information(this, "!!", QString::number(nChart));
        chartLayout->addWidget(buttons[ii+stChart], ii+1, 0);
        chartLayout->addWidget(chartview[ii+stChart], ii+1, 1);
    }
    this->setLayout(chartLayout);

}

void QWidgetWAVChart::setData(){
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    QVector<double>& wavt = parentwindow->getDatat();
    for (int ii=0; ii<nChart; ++ii){
        axesX[ii]->setRange(wavt[0],wavt[wavt.length()-1]);
        //QMessageBox::information(this, "!!", QString::number(electrodearray[ii+stChart]));
        //QMessageBox::information(this, "!!", QString::number(electrodearray[ii+stChart]) + ":" + channelname[electrodearray[ii+stChart]]);
        axesX[ii]->setTitleText('[' + QString::number(ii+1+stChart)  + ']' + "wavraw");
    }
    updateSeriesX();
    saveAxesXHistory();

    //QMessageBox::information(this, "!!", "ok");
    calcSpectrogram();
}

void QWidgetWAVChart::calcSpectrogram(){
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    QVector<qint16>& wavraw = parentwindow->getWavraw();
    QVector<double>& wavt = parentwindow->getDatat();

    int freq = 44100;
    double winwid = 0.01;
    double padwid = 0.01;
    double shiftwid = 0.001;
    double stt = 0;
    double ent = wavt[wavt.length()-1];
    int stpnt = int(stt*freq); // starting point of calculation
    int enpnt = int(ent*freq-1); // end point of calculation
    if (stpnt < 0){
        stpnt = 0;
    }
    if (enpnt > wavt.length()-1){
        enpnt = wavt.length()-1;
    }

    int winpntorg = int(winwid*44100); // winwid*44100
    int padpnt = int(padwid*44100);
    int shiftpnt = int(shiftwid*44100); // this truncation will change the calc times

    int max_freq = 1*44100;
    int min_freq = 1/winwid;
    double padCoef = padwid/winwid;
    double dpnt = winwid/shiftwid;
    double df;
    if (padCoef > 1){
        df = min_freq/padCoef;
    }else{
        df = min_freq;
    }

    // These values have to be even value
    if (winpntorg%2 != 0) winpntorg++;
    if (padpnt%2 != 0) padpnt++;
    int winpnt = padpnt;

    int Ncalc = int((enpnt - stpnt - winpnt)/shiftpnt)+1;
    Spec_mag.resize(Ncalc);
    Spec_phase.resize(Ncalc);
    Spec_x.resize(Ncalc);

    //QMessageBox::information(this, "Unable to save file", QString::number(int((enpnt - stpnt - winpnt)/shiftpnt))+1);

    // Use FFTW library
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * winpnt);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * winpnt);

    // Use this if you do FFT for a fixed N data multiple times (it takes a few sec for init but then it uses the fastest algorithm)
    fftw_plan p = fftw_plan_dft_1d(winpnt, in, out, FFTW_FORWARD, FFTW_MEASURE);
    // Or, use this if you do FFT for multiple data with different size (it doesn't need init time, but may not be optimal)
    //fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    Spec_mag_max = 0;

    double hanning = 1;

    int sindex = 0;
    int stwin = stpnt;
    int enwin = stwin + winpnt;
    do {
        // Prepare data
        for(int ii=0; ii<winpnt; ii++){
            if (ii < (padpnt-winpnt)/2 || (padpnt+winpnt)/2 <= ii){
                in[ii][0] = 0; // real value, padding
            } else {
                hanning = 0.5 * (1 - cos(2 * PI * (ii-(padpnt-winpnt)/2) / (winpnt+(padpnt-winpnt)/2-1)));
                in[ii][0] = double(wavraw[stwin+ii]) * hanning; // real value
            }
            in[ii][1] = 0; // imaginary value
        }
        fftw_execute(p); // FFT (real output: out[ii][0]; imaginary output: out[ii][1])

        // Spec save
        Spec_mag[sindex].resize(winpnt/2); // magnitude in [quantity peak]
        Spec_phase[sindex].resize(winpnt/2); // phase in [radians] (-pi ~ pi)
        for(int ii=0; ii<winpnt/2; ii++){
            if (ii == 0 || ii == winpnt/2){ // or you can include the Nyquist (Nyquist is discarded in the current code)
                Spec_mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt; // DC and Nyquist ([rms] magnitude is the same)
                //FFTmaglog[ii] = 20 * log(10,FFTmag[ii])/; // powerspectrum
                //FFTpower[ii] = FFTmag[ii]*FFTmag[ii]; // powerspectrum
            } else {
                Spec_mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt * 2; // pos+neg frequencies (further dividing with sqrt(2) will give you [rms] magnitude)
                //FFTpower[ii] = FFTmag[ii]*FFTmag[ii]/2; // powerspectrum
            }
            if (Spec_mag_max < Spec_mag[sindex][ii]){
                Spec_mag_max = Spec_mag[sindex][ii];
            }
            Spec_phase[sindex][ii] = atan(out[ii][1] / out[ii][0]); // DC and Nyquist
        }
        Spec_x[sindex] = stt + winwid/2 + double(sindex * shiftpnt) / freq;

        sindex++;
        stwin = stpnt + sindex * shiftpnt;
        enwin = stwin + winpnt;
    } while(sindex < Ncalc-1);
    //QMessageBox::information(this, "Unable to save file", QString::number(sindex));

    Spec_y.resize(winpnt/2);
    for(int ii=0; ii<winpnt/2; ii++){
        Spec_y[ii] = 0 + ii*df;
    }
    Spec_dx = shiftwid;
    Spec_dy = df;

    // Destroy FFTW
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    prepareSpecImage();
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
    for (int ii=0; ii<nChart; ++ii){
        series[ii]->clear();
    }

    QVector<qint16>& wavraw = parentwindow->getWavraw();
    QVector<double>& wavt = parentwindow->getDatat();

    for (int ii=0; ii<nChart; ++ii){
    // for (int ii=stChart; ii<nChart+stChart; ++ii){
    // *** ii = #Ch
        qreal minaxisX = axesX[ii]->min();
        qreal maxaxisX = axesX[ii]->max();

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
        //QMessageBox::information(this, "!!tlength", "Ch:" + QString::number(electrodearray[ii+stChart]));

        minx = wavt[stpnt];
        miny = wavraw[stpnt];
        maxx = wavt[stpnt];
        maxy = wavraw[stpnt];
        int allminy = int(miny);
        int allmaxy = int(maxy);

        if (wavtlen <= nPoints){
            for (int jj = stpnt; jj < enpnt; ++jj){
                series[ii]->append(wavt[jj], wavraw[jj]);
                //seriesfilt[ii]->append(wavt[jj], datayfilt[electrodearray[ii]][jj]);
                if (allminy > wavraw[jj]){
                    allminy = int(wavraw[jj]);
                }
                if (allmaxy < wavraw[jj]){
                    allmaxy = int(wavraw[jj]);
                }
            }
        } else {
            for (int jj = stpnt; jj < enpnt; ++jj){
                if (jj%(wavtlen/nPoints)==0){
                    if (minx < maxx){
                        series[ii]->append(minx, miny);
                        series[ii]->append(maxx, maxy);
                        //seriesfilt[ii]->append(minx, minyfilt);
                        //seriesfilt[ii]->append(maxx, maxyfilt);
                    } else {
                        series[ii]->append(maxx, maxy);
                        series[ii]->append(minx, miny);
                        //seriesfilt[ii]->append(maxx, maxyfilt);
                        //seriesfilt[ii]->append(minx, minyfilt);
                    }
                    if (allminy > miny){
                        allminy = (int)miny;
                    }
                    if (allmaxy < maxy){
                        allmaxy = (int)maxy;
                    }
                    miny = wavraw[jj];
                    maxy = wavraw[jj];
                    minx = wavt[jj];
                    maxx = wavt[jj];
                    //maxyfilt = datayfilt[electrodearray[ii]][jj];
                    //minyfilt = datayfilt[electrodearray[ii]][jj];
                } else {
                    if (maxy < wavraw[jj]){
                        maxy = wavraw[jj];
                        maxx = wavt[jj];
                        //maxyfilt = datayfilt[electrodearray[ii]][jj];
                    }
                    if (miny > wavraw[jj]){
                        miny = wavraw[jj];
                        minx = wavt[jj];
                        //minyfilt = datayfilt[electrodearray[ii]][jj];
                    }
                }
            }

        }

        if (axesYmin[ii] == 0 && axesYmax[ii] == 0){
            axesY[ii]->setRange((allminy - allminy%100 - 100), (allmaxy - allmaxy%100 + 100));
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


    //QMessageBox::information(this, "!!tlength", QString::number(wavt.length()));
    //QMessageBox::information(this, "!!stpnt", QString::number(stpnt));
    //QMessageBox::information(this, "!!enpnt", QString::number(enpnt));

    prepareSpecImage();
    repaint();
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
    for (int ii=0; ii<nChart; ++ii){
        axesX[ii]->setRange(min,max);
    }
    updateSeriesX();
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
    for (int ii=0; ii<nChart; ++ii){
        qreal min = axesX[ii]->min();
        qreal max = axesX[ii]->max();
        qreal len = max - min;
        axesX[ii]->setRange(min + len*scalemin, max + len*scalemax);
    }
    updateSeriesX();
}

void QWidgetWAVChart::alignAxesXto(int indChart){
    qreal min = axesX[indChart]->min();
    qreal max = axesX[indChart]->max();
    for (int ii=0; ii<nChart; ++ii){
        axesX[ii]->setRange(min,max);
    }
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
        historymin[9] = axesX[0]->min();
        historymax[9] = axesX[0]->max();
    } else {
        if (historyind == historymin.size()){
            historymin.append(axesX[0]->min());
            historymax.append(axesX[0]->max());
            historyind += 1;
        } else {
            historymin[historyind] = axesX[0]->min();
            historymax[historyind] = axesX[0]->max();
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
    if (Spec_x.length() == 0) return;
    //QMessageBox::information(this, "!prepareSpecImage", "ok");
    axesY[1]->setRange(0, 10000);
    // left & right t
    qreal minx = axesX[1]->min(); // mint ~ Topleft.x()
    qreal maxx = axesX[1]->max(); // maxt ~ Topleft.x() + width
    qreal miny = axesY[1]->min(); // mint ~ Topleft.x()
    qreal maxy = axesY[1]->max(); // maxt ~ Topleft.x() + width
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
    specoffsety /= xwid;
    specboxx = Spec_dx/xwid;
    specboxy = Spec_dy/ywid;

    //QMessageBox::information(this, "!specoffsetx", QString::number(specoffsetx));
    //QMessageBox::information(this, "!specoffsety", QString::number(specoffsety));
    //QMessageBox::information(this, "!specboxx", QString::number(specboxx));
    //QMessageBox::information(this, "!specboxy", QString::number(specboxy));
}



void QWidgetWAVChart::paintEvent(QPaintEvent* event){
    if (specnumx == 0 || specnumy == 0) return;

    // prepare image space for spectrogram
    qreal plotwidth = chart[0]->plotArea().width();
    qreal plotheight = chart[0]->plotArea().height();
    qreal viewwidth = chartview[0]->width();
    qreal viewheight = chartview[0]->height();    // Spectrograms
    QImage imgspec_bg(viewwidth, viewheight, QImage::QImage::Format_ARGB32);
    imgspec_bg.fill(Qt::white); // assuming this is white
    QImage imgspec(plotwidth, plotheight, QImage::Format_ARGB32); // check pixel manipulation
    QPainter painter(&imgspec_bg);
    QPointF TopLeft = chart[1]->plotArea().topLeft();
    painter.drawImage(TopLeft, imgspec);

    // paint spectrogram (this is not optimal because the calc is more than resolution)
    imgspec.fill(Qt::white);
    QPainter mypainter(&imgspec);
    QRect tmprec;
    int gainshade = 50;
    int colorval = 0;
    for (int ii=0; ii<specnumx; ii++){
        tmprec.setLeft(TopLeft.x() + int((specoffsetx - specboxx/2 + specboxx * ii)*plotwidth));
        tmprec.setWidth(int(specboxx * plotwidth)+1);
        for (int jj=0; jj<specnumy; jj++){
            tmprec.setTop(TopLeft.y() + int((1 - (specoffsety + specboxy/2 + specboxy * jj))*plotheight));
            tmprec.setHeight(int(specboxy * plotheight)+1);

            colorval = 255 - int(Spec_mag[ii+specstindx][jj+specstindy]/Spec_mag_max * 255) * gainshade;
            if (colorval < 0) colorval = 0;
            mypainter.fillRect(tmprec, QColor(colorval, colorval, colorval));
//            mypainter.fillRect(tmprec, Qt::red);
        }
    }
    //tmprec.setLeft(int(TopLeft.x()));
    //tmprec.setTop(int(TopLeft.y()));
    //tmprec.setWidth(int(0.002 * plotwidth));
    //tmprec.setHeight(int(0.002 * plotheight));
    //mypainter.fillRect(tmprec, Qt::red);
    chart[1]->setPlotAreaBackgroundBrush(imgspec);
    chart[1]->setPlotAreaBackgroundVisible(true);


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
