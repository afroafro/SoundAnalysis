#include "mainwindow.h"
#include "qwidgetwavchart.h"
#include "ui_mainwindow.h"
#include <spuc/generic/iir_2nd.h>
#include <spuce/filters/iir_coeff.h>
#include <spuce/filters/elliptic_iir.h>
#include <thread>
#include "../../fftw-3.3.5-dll32/fftw3.h"

using namespace SPUC;
using namespace spuce; // removed several files: filter::audio_equalizer.cpp/h, create_remez_fir.cpp

/*!
    \class MainWindow
            This class generates mainwindow and related ui components, which can open WAV files and analyze or save the data.
*/

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    // parameters (these may not be working; look at WidgetIntanChart.cpp)
    fileopened = false;
    appbusy = false;
    /*/////////////////////////////////////*/

    specs.resize(1);    // without this, qwidgetwavchart generates error during repaint()

    // *** Create widget
    widget = new QWidget(this);
    mainLayout = new QGridLayout;
    QMenuBar* menuBar = new QMenuBar;
    QMenu* fileMenu = new QMenu("&File", widget);
    QMenu* displayMenu = new QMenu("&Display", widget);
    menuBar->addMenu(fileMenu);
    menuBar->addMenu(displayMenu);
    const int nfileaction = 10;
    QAction* fileAction[nfileaction];
    fileAction[0] = new QAction("&Open wav");
    fileAction[1] = new QAction("Open2");
    fileAction[2] = new QAction("&Save Data");
    fileAction[3] = new QAction("Save TutorData");
    fileAction[4] = new QAction("Save spike time");
    fileAction[5] = new QAction("Load Data");
    fileAction[6] = new QAction("Load TutorData");
    fileAction[7] = new QAction("Auto sort");
    fileAction[8] = new QAction("Auto extraction");
    fileAction[9] = new QAction("Delete TutorData");
    for (int ii=0; ii<nfileaction; ++ii){
        fileAction[ii]->setData(ii);
        fileMenu->addAction(fileAction[ii]);
    }
    const int ndisplayaction = 4;
    QAction* displayAction[ndisplayaction];
    displayAction[0] = new QAction("Reorder Channels");
    displayAction[1] = new QAction("Show consecutive files");
    displayAction[2] = new QAction("Show variables");
    displayAction[3] = new QAction("Show sorted spikes");
    for (int ii=0; ii<ndisplayaction; ++ii){
        displayAction[ii]->setData(ii);
        displayMenu->addAction(displayAction[ii]);
    }
    mainLayout->setMenuBar(menuBar);
    connect(fileMenu, SIGNAL(triggered(QAction*)), this, SLOT(pushFile(QAction*)));
    connect(displayMenu, SIGNAL(triggered(QAction*)), this, SLOT(pushDisplay(QAction*)));
    widget->setLayout(mainLayout);
    this->setCentralWidget(widget);

    widgetchart0 = new QWidgetWAVChart(this);
    mainLayout->addWidget(widgetchart0, 0, 0);
    bottomwidget = new QWidget(this);
    mainLayout->addWidget(bottomwidget, 1, 0);

    // bottom
    bottomLayout = new QGridLayout;
    //mainLayout->addLayout(bottomLayout,1,0);
    bottomwidget->setLayout(bottomLayout);
    listwidget = new QListWidget(this);
    //listwidget->setMinimumWidth(listwidget->sizeHintForColumn(0));
    listwidget->addItem("Open any directory");
    bottomLayout->addWidget(listwidget, 0, 0);
    buttonwidget = new QWidget(this);
    bottomLayout->addWidget(buttonwidget, 0, 1);

    QObject::connect(listwidget, SIGNAL(itemSelectionChanged()), this, SLOT(openFile()));

    //bottomLayout->addWidget(buttonwidget, 0, 1);
    //bottomLayout->addWidget(listwidget1, 0, 1);

    // *** initialize widget
    //widgetchart0->initParams();

    // Buttons
    buttons.append(new QPushButton("Amp"));
    buttons.append(new QPushButton("Filter"));
    buttons.append(new QPushButton("Waves"));
    buttons.append(new QPushButton("Batch detection"));
    buttons.append(new QPushButton("Stop batch"));
    buttons.append(new QPushButton("Detect sound"));
    buttons.append(new QPushButton("Play sound"));
    buttons.append(new QPushButton("Stop sound"));
    buttons.append(new QPushButton("Play sound part"));
    buttons.append(new QPushButton("Save WAV"));
    buttons.append(new QPushButton("Rhythm ana"));
    QObject::connect(buttons[0], SIGNAL(released()), this, SLOT(pushAmp()));
    QObject::connect(buttons[1], SIGNAL(released()), this, SLOT(pushFilter()));
    QObject::connect(buttons[2], SIGNAL(released()), this, SLOT(pushWaves()));
    QObject::connect(buttons[3], SIGNAL(released()), this, SLOT(pushBatch()));
    QObject::connect(buttons[4], SIGNAL(released()), this, SLOT(pushStopBatch()));
    QObject::connect(buttons[5], SIGNAL(released()), this, SLOT(detectSound()));
    QObject::connect(buttons[6], SIGNAL(released()), this, SLOT(pushPlaySound()));
    QObject::connect(buttons[7], SIGNAL(released()), this, SLOT(pushStopSound()));
    QObject::connect(buttons[8], SIGNAL(released()), this, SLOT(pushPlaySoundPart()));
    QObject::connect(buttons[9], SIGNAL(released()), this, SLOT(pushSaveWAV()));
    QObject::connect(buttons[10], SIGNAL(released()), this, SLOT(pushRhythm()));

    QGridLayout* buttonLayout = new QGridLayout;
    buttonLayout->setAlignment(Qt::AlignTop);
    buttonwidget->setLayout(buttonLayout);
    for (int ii=0; ii<buttons.length(); ii++){
        buttonLayout->addWidget(buttons[ii], ii, 0);
    }

    specslider = new QSlider(Qt::Vertical);
    specslider->setFocusPolicy(Qt::StrongFocus);
    specslider->setTickPosition(QSlider::TicksBothSides);
    specslider->setTickInterval(50);
    specslider->setSingleStep(1);
    specslider->setMinimum(1);
    specslider->setMaximum(100);
    specslider->setValue(50);
    connect(specslider, SIGNAL(valueChanged(int)), this, SLOT(specsliderChanged()));
    buttonLayout->addWidget(specslider, 0, 1);

    //player = new QMediaPlayer(this, QMediaPlayer::StreamPlayback);

    disableButtons();
}

MainWindow::~MainWindow() {

    /*
    if (buffer->isOpen()){ // this causes crash
        buffer->close();
    }
    */
    delete ui;
}

void MainWindow::pushFile(QAction *action){
    int value = action->data().toInt();
    if (value == 0){
        openDir();
    } else if (value == 1){
    } else if (value == 2){
        widgetchart0->saveData();
    } else if (value == 3){
    } else if (value == 4){
    } else if (value == 5){
    } else if (value == 6){
    } else if (value == 7){
    } else if (value == 8){
    } else if (value == 9){
    }
}

void MainWindow::openFile(){
    // Opens .wav file and initialize all the data including spectrograms
    // Some notification before initializing is desirable to prevent data loss

    // This might be better: "https://www.ipentec.com/document/csharp-play-wave-file-using-wave-api"
    if (appbusy || wav.filePath.isEmpty()){
        return;
    }

    QString fileNametmp = listwidget->currentItem()->text();
    QString filePathtmp = wav.filePath + "/" + fileNametmp;

    //QMessageBox::information(this, "file path", filePathtmp);
    //QMessageBox::information(this, "file name", fileNametmp);

    if (fileNametmp.isEmpty()){
        return;
    } else {
        appbusy = true;
        this->setCursor(Qt::WaitCursor);
        QFile file(filePathtmp);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(this, "Unable to open file", file.errorString());
            appbusy = false;
            this->setCursor(Qt::ArrowCursor);
            return;
        }
        QDataStream in(&file);
        in.setVersion(QDataStream::Qt_4_8); // should I use pointer? in -> setVersion ?
        in.setByteOrder(QDataStream::LittleEndian);
        in.setFloatingPointPrecision(QDataStream::SinglePrecision);

        wav.filebytesize = file.size();

        in >> wav.RIFF[0] >> wav.RIFF[1] >> wav.RIFF[2] >> wav.RIFF[3];
        in >> wav.wavfilesize;
        in >> wav.WAVE[0] >> wav.WAVE[1] >> wav.WAVE[2] >> wav.WAVE[3];
        in >> wav.fmt[0] >> wav.fmt[1] >> wav.fmt[2] >> wav.fmt[3];
        in >> wav.ChunkSize; // 16
        in >> wav.AudioFormat ;
        in >> wav.NCh;
        in >> wav.SamplesPerSec;
        in >> wav.bytesPerSec;
        in >> wav.blockAlign ;
        in >> wav.bitsPerSample;
        in >> wav.Subchunk2ID[0] >> wav.Subchunk2ID[1] >> wav.Subchunk2ID[2] >> wav.Subchunk2ID[3];
        in >> wav.Subchunk2Size;
        wav.datasamplesize = wav.Subchunk2Size / (wav.NCh * wav.bitsPerSample/8);
        if (wav.NCh == 1){
            wav.dataL.resize(wav.datasamplesize);
            for(int ii=0; ii<wav.dataL.size(); ii++){
                in >> wav.dataL[ii];
            }
        } else if (wav.NCh == 2){
            wav.dataL.resize(wav.datasamplesize);
            wav.dataR.resize(wav.datasamplesize);
            for(int ii=0; ii<wav.dataL.size(); ii++){
                in >> wav.dataL[ii];
                in >> wav.dataR[ii];
            }
        } else {
            // error
            file.close();
            QMessageBox::information(this, "!Exception", "Seems not match any data format.");
            appbusy = false;
            this->setCursor(Qt::ArrowCursor);
            return;
        }
        file.close();

        QString msg = "";
        msg += "RIFF:\t" + QString::number(wav.RIFF[0]) + QString::number(wav.RIFF[1]) + QString::number(wav.RIFF[2]) + QString::number(wav.RIFF[3]) + "\n";
        msg += "wavfile:\t" + QString::number(wav.wavfilesize) + "\n";
        msg += "WAVE:\t" + QString::number(wav.WAVE[0]) + QString::number(wav.WAVE[1]) + QString::number(wav.WAVE[2]) + QString::number(wav.WAVE[3]) + "\n";
        msg += "fmt:\t" + QString::number(wav.fmt[0]) + QString::number(wav.fmt[1]) + QString::number(wav.fmt[2]) + QString::number(wav.fmt[3]) + "\n";
        msg += "Chunksize:\t" + QString::number(wav.ChunkSize) + "\n";
        msg += "audio:\t" + QString::number(wav.AudioFormat) + "\n";
        msg += "NCh:\t" + QString::number(wav.NCh) + "\n";
        msg += "samples/s:\t" + QString::number(wav.SamplesPerSec) + "\n";
        msg += "bytes/s:\t" + QString::number(wav.bytesPerSec) + "\n";
        msg += "blockalign:\t" + QString::number(wav.blockAlign) + "\n";
        msg += "bits/sample:\t" + QString::number(wav.bitsPerSample) + "\n";
        msg += "subchunk:\t" + QString::number(wav.Subchunk2ID[0]) + QString::number(wav.Subchunk2ID[1]) + QString::number(wav.Subchunk2ID[2]) + QString::number(wav.Subchunk2ID[3]) + "\n";
        msg += "subchunksize:\t" + QString::number(wav.Subchunk2Size) + "\n";
        //QMessageBox::information(this, "!Done", msg);
    }

    wav.fileName = fileNametmp;
    //wav.filePath = filePathtmp;
    QStringList dirnames = filePathtmp.split("/");
    wav.fileDir = "";
    for (int ii=0; ii<dirnames.length()-1; ++ii){
        wav.fileDir = wav.fileDir + dirnames[ii] + "/";
    }

    fileopened = true;

    wav.wavt.resize(wav.dataL.length());
    for (int ii=0; ii<wav.dataL.length(); ii++){
        wav.wavt[ii] = double(ii)/double(wav.SamplesPerSec);
    }

    // Making v[0] as the copy of wav.dataL (but is it possible that dataL is empty??? then we will need to throw some exceptions)
    if (wav.dataL.length() > 0){
        v.resize(1);
        v[0].data.resize(wav.dataL.size());
        for(int ii=0; ii<wav.dataL.size(); ii++){
            v[0].data[ii] = double(wav.dataL[ii]);
        }
        v[0].samplingfreq = wav.SamplesPerSec;
    }
    if (wav.dataR.length() > 0){ // R data is stored in another vector
        v.resize(v.size()+1);
        v[v.size()-1].data.resize(wav.dataR.size());
        for(int ii=0; ii<wav.dataR.size(); ii++){
            v[v.size()-1].data[ii] = double(wav.dataR[ii]);
        }
        v[v.size()-1].samplingfreq = wav.SamplesPerSec;
    }
    init();

    //player.setMedia(QUrl::fromLocalFile(filePathtmp));
    //player.setPlaybackRate(1.0);
    //player.setVolume(100);

    this->setWindowTitle("masaWavRead: " + wav.fileName);
    appbusy = false;
    this->setCursor(Qt::ArrowCursor);

    widgetchart0->setData();

    specs[0].ind = 0;
    specs[0].winwid = 0.01;
    specs[0].padwid = 0.01;
    specs[0].shiftwid = 0.005;
    widgetchart0->setAxisYspec(0, 0, 20000);
    calcSpectrogram(v[0].data, wav.wavt, specs[0]);

    widgetchart0->prepareSpecImage();

    enableButtons();
}
void MainWindow::enableButtons(){
    for (int ii=0; ii<buttons.length(); ii++){
        buttons[ii]->setEnabled(true);
    }
}

void MainWindow::disableButtons(){
    for (int ii=0; ii<buttons.length(); ii++){
        buttons[ii]->setEnabled(false);
    }
}


void MainWindow::pushAmp(){
    for (int ii=0; ii<wav.dataL.length(); ii++){
        wav.dataL[ii] = 0;
    }
    //widgetchart0->setSeries();
}

void MainWindow::pushFilter(){
    QVector<QString> dialoglabel;
    QVector<QString> dialogvalue;
    QVector<QString> returnvalue;
    // other potential params: stopband?
    dialoglabel.append("Highpass cut-off freq [Hz] (skipped if<0)");
    dialoglabel.append("Highpass attenuation [dB]");
    dialoglabel.append("Highpass ripples [dB]");
    dialoglabel.append("Lowpass cut-off freq [Hz] (skipped if<0)");
    dialoglabel.append("Lowpass attenuation [dB]");
    dialoglabel.append("Lowpass ripples [dB]");
    dialogvalue.append(QString::number(highpass_freq));
    dialogvalue.append(QString::number(highpass_attn));
    dialogvalue.append(QString::number(highpass_ripple));
    dialogvalue.append(QString::number(lowpass_freq));
    dialogvalue.append(QString::number(lowpass_attn));
    dialogvalue.append(QString::number(lowpass_ripple));
    returnvalue = myDialog("Elliptic filter (2nd-order):", dialoglabel, dialogvalue);
    if (returnvalue.length() == 0){
        return;
    }
    highpass_freq = returnvalue[0].toDouble();
    highpass_attn = returnvalue[1].toDouble();
    highpass_ripple = returnvalue[2].toDouble();
    lowpass_freq = returnvalue[3].toDouble();
    lowpass_attn = returnvalue[4].toDouble();
    lowpass_ripple = returnvalue[5].toDouble();

    struct FiltParam fp;
    fp.samplerate = wav.SamplesPerSec;
    fp.HPfreq = highpass_freq;
    fp.HPattn = highpass_attn;
    fp.HPripple = highpass_ripple;
    fp.LPfreq = lowpass_freq;
    fp.LPattn = lowpass_attn;
    fp.LPripple = lowpass_ripple;

    QVector<double> wave0;
    wave0 = filterWave(v[0].data, fp);
    addAnaWave(wave0);
}

void MainWindow::addAnaWave(QVector<double> wave0){
    v.resize(v.size()+1);
    v[v.size()-1].data = wave0;
    int anaind = v.size()-1;
    widgetchart0->setSeries(anaind, 1);
}


QVector<double> MainWindow::filterWave (QVector<double> wave0, struct FiltParam fp){

    uint iir_coeff_order = 2;
    if (fp.HPfreq <= 0) fp.HPfreq = 0.001; // Hz
    if (fp.LPfreq <= 0) fp.LPfreq = 0.001;
    double highpassfreq = fp.HPfreq/double(fp.samplerate); // cut-off freq (1=sampling rate)
    double lowpassfreq = fp.LPfreq/double(fp.samplerate); // cut-off freq (1=sampling rate)
    double highpassattn = fp.HPattn; // passband attenuation. dB?
    double lowpassattn = fp.LPattn; // passband attenuation. dB?
    double highpassripple = fp.HPripple; // ripple. dB?
    double lowpassripple = fp.LPripple; // ripple. dB?

    //highpass filter
    iir_coeff* filt = new iir_coeff(iir_coeff_order, filter_type::high);
    if (highpassfreq <= 0){
        elliptic_iir(*filt, 0.000000001, 0.000000001, 0.000000001);
    } else if (highpassfreq >= 0.5){
        elliptic_iir(*filt, 0.499999999, highpassripple, highpassattn);
    } else {
        elliptic_iir(*filt, highpassfreq, highpassripple, highpassattn);
    }
    iir_2nd<double> highpassfilt(filt->get_b(0), filt->get_b(1), filt->get_b(2), filt->get_a(1), filt->get_a(2));

    // lowpass filter
    filt = new iir_coeff(iir_coeff_order, filter_type::low);
    if (lowpassfreq <= 0){
        elliptic_iir(*filt, 2, 0.000000001, 0.000000001);
    } else if (lowpassfreq >= 0.5){
        elliptic_iir(*filt, 0.499999999, lowpassripple, lowpassattn);
    } else {
        elliptic_iir(*filt, lowpassfreq, lowpassripple, lowpassattn);
    }
    iir_2nd<double> lowpassfilt(filt->get_b(0), filt->get_b(1), filt->get_b(2), filt->get_a(1), filt->get_a(2));

    QVector<double> tmpwave;
    tmpwave.resize(wave0.length());
    for (int ii=0; ii<wave0.length(); ++ii){
        tmpwave[ii] = highpassfilt.clock(lowpassfilt.clock(wave0[ii]));
    }
    delete filt;

    return tmpwave;
}

QVector<double> MainWindow::absWave(QVector<double> wave0){
    QVector<double> tmpwave;
    tmpwave.resize(wave0.length());
    for (int ii=0; ii<wave0.length(); ++ii){
        if (wave0[ii]<0){
            tmpwave[ii] = -wave0[ii];
        } else {
            tmpwave[ii] = wave0[ii];
        }
    }
    return tmpwave;
}

QVector<double> MainWindow::differentiateWave(QVector<double> wave0){
    QVector<double> tmpwave;
    tmpwave.resize(wave0.length()-1);
    for (int ii=0; ii<tmpwave.length(); ++ii){
        tmpwave[ii] = wave0[ii+1]-wave0[ii];
    }
    return tmpwave;
}

QVector<double> MainWindow::avgDiscretizeWave(QVector<double> wave0, int step){
    int steplen = int(wave0.length()/step);

    QVector<double> tmpwave;
    tmpwave.resize(steplen);
    for (int ii=0; ii<steplen; ++ii){
        double sum = 0.0;
        for (int jj=0; jj<step; ++jj){
            sum += wave0[ii*step+jj];
        }
        tmpwave[ii] = sum / step;
    }
    return tmpwave;
}

void MainWindow::pushDisplay(QAction *action){
    int value = action->data().toInt();
    if (value == 0){
    } else if (value == 1){
    } else if (value == 2){ // show parameters
    } else if (value == 3){ // show sorted spikes
    }
}

void MainWindow::pushWaves(){
    QVector<QVector<int>> axesYseries;
    axesYseries = widgetchart0->getAxesYseries();
    QVector<int> seriesYaxes;
    seriesYaxes.resize(v.length());
    for (int kk=0; kk<v.length(); kk++){
        for (int ii=0; ii<axesYseries.length(); ii++){
            for (int jj=0; jj<axesYseries[ii].length(); jj++){
                if (kk == axesYseries[ii][jj]){
                    seriesYaxes[kk] = ii;
                }
            }
        }
    }
    QString msg = "";
    for (int kk=0; kk<v.length(); kk++){
        msg += "Series[" + QString::number(kk) + "]: axisY" + QString::number(seriesYaxes[kk]) + "\n";
    }
    QMessageBox::information(this, "!Done", msg);
}

void MainWindow::pushRhythm(){
    v.resize(v.size()+1);
    int anaind = v.size()-1;
    v[anaind].data.resize(wav.dataL.length());

    // Band pass filter 200-20000Hz
    struct FiltParam fp;
    fp.samplerate = wav.SamplesPerSec;
    fp.HPfreq = 200;
    fp.HPattn = 50;
    fp.HPripple = 0.1;
    fp.LPfreq = 20000;
    fp.LPattn = 50;
    fp.LPripple = 0.1;
    v[anaind].data = filterWave(v[0].data, fp);

    // Rectification
    v[anaind].data = absWave(v[anaind].data);

    // Discretize with moving average of step 10pnt (10/44100 s) 0.01 * 44100 = 441
    double stept = 0.01; // [s]
    int steppnt = int(stept * wav.SamplesPerSec);
    v[anaind].data = avgDiscretizeWave(v[anaind].data, steppnt);
    QVector<double> wavtdiscretize;
    wavtdiscretize = avgDiscretizeWave(wav.wavt, steppnt);

    // log-transform wave
    for (int ii=0; ii<v[anaind].data.length(); ++ii){
        v[anaind].data[ii] = log(v[anaind].data[ii])/log(10)*10;
    }

    // differentiate wave
    for (int ii=0; ii<v[anaind].data.length()-1; ++ii){
        v[anaind].data[ii] = v[anaind].data[ii+1]-v[anaind].data[ii];
        wavtdiscretize[ii] = (wavtdiscretize[ii+1]+wavtdiscretize[ii])/2;
    }
    v[anaind].data.resize(v[anaind].data.length()-1);
    wavtdiscretize.resize(wavtdiscretize.length()-1);
    v[anaind].samplingfreq = round(double(wav.SamplesPerSec)/double(steppnt));

    // LP filter < 20Hz
    struct FiltParam fp2;
    fp2.samplerate = round(double(wav.SamplesPerSec)/double(steppnt));
    fp2.HPfreq = 0;
    fp2.HPattn = 50;
    fp2.HPripple = 0.1;
    fp2.LPfreq = 20;
    fp2.LPattn = 50;
    fp2.LPripple = 0.1;
    v[anaind].data = filterWave(v[anaind].data, fp2);
    //saveData1D(v[anaind].data);

    specs[0].ind = anaind;
    specs[0].winwid = 1;
    specs[0].padwid = 4;
    specs[0].shiftwid = 0.02;
    calcSpectrogram(v[anaind].data, wavtdiscretize, specs[0]); //t=0.012483 (0) - 15886.5 (1588648)
    widgetchart0->setAxisYspec(0, 0, 50);
    widgetchart0->prepareSpecImage();

    //saveData2D(specs[0].mag);
}





void MainWindow::pushBatch(){
    int indnow = listwidget->currentRow();
    userstopped = false;
    while (indnow+1 < listwidget->count() && !userstopped){
        listwidget->setCurrentRow(indnow+1);
        openFile();
        detectSound();
        indnow++;
        QCoreApplication::processEvents();
    }
}

void MainWindow::pushStopBatch(){
    userstopped = true;
}

void MainWindow::detectSound(){
    pushFilter();
    double winsize = 0.005;
    double winstepsize = 1.0 / wav.SamplesPerSec;
    double sampling = wav.SamplesPerSec;

    double winstep = round(winstepsize * sampling);
    double avgbase=0.0, sdbase=0.0, varbase=0.0;
    double pntwinsize = round(winsize * sampling);
    bool flagbase = true;
    int nbase = 0;
    int lenwav = v[1].data.length();
    double thres = 3.0;
    double minDur = 20.0; //[ms]
    double minISI = 1.0; // [ms]

    double avg, sd, M2, tmpavg, delta;
    for(int ii=0; ii<lenwav-pntwinsize+1; ii+=winstep){
        avg = 0.0;
        sd = 0.0;
        delta = 0.0;
        M2 = 0.0;
        for(int jj=0; jj<pntwinsize; jj++){ // this algorithm is numerically instable. use Welford's algorithm instead
            //delta = anaDataL[0][ii+jj]-avg;
            //avg += delta/pntwinsize;
            //M2 += delta*(anaDataL[0][ii+jj]-avg);

            tmpavg = avg;
            avg += v[1].data[ii+jj];
            sd += v[1].data[ii+jj] * v[1].data[ii+jj];
            //sd += ((anaDataL[0][ii+jj]-tmpavg)*(anaDataL[0][ii+jj]-avg)-sd)/pntwinsize;
            //QMessageBox::information(this, "file path", "anadata: [" + QString::number(jj) + "]" + QString::number(anaDataL[0][ii+jj]));
        }
        avg /= pntwinsize;
        sd = sqrt((sd/pntwinsize - avg*avg)*pntwinsize/(pntwinsize-1));
        //sd = sqrt(M2/(pntwinsize-1));
        //QMessageBox::information(this, "file path", "st: " + QString::number(ii) + "\r\n en: " + QString::number(ii+pntwinsize) + "\r\n avg: " + QString::number(avg) + "\r\n sd: " + QString::number(sd));

        if (ii==0){
            avgbase = sd;
            varbase = sd*sd;
            nbase += 1;
            continue;
        }
        if (ii==1){
            avgbase = avgbase*nbase/(nbase+1) + sd/(nbase+1);
            //varbase = varbase*nbase + V_sdev*V_sdev // twice larger
            varbase = varbase*nbase/(nbase+1) + sd*sd/(nbase+1);
            sdbase = sqrt(varbase - avgbase*avgbase);
            nbase += 1;
            continue;
        }

        //QMessageBox::information(this, "file path", QString::number(ii) + "\r\n" + QString::number(sd) + "\r\n" + QString::number(avgbase) + "\r\n" + QString::number(sdbase));
        //QMessageBox::information(this, "file path", QString::number(pntwinsize) + "\r\n" + QString::number(ii) + "\r\n" + QString::number(sd) + "\r\n" + QString::number(avgbase) + "\r\n" + QString::number(sdbase));

        if (sd > avgbase + thres * sdbase){
            if (flagbase){
                flagbase = false;
                ton.append(wav.wavt[ii+int(pntwinsize/2)]);
            }
        } else {
            if (!flagbase){
                flagbase = true;
                toff.append(wav.wavt[ii+int(pntwinsize/2)]);
            }
            avgbase = avgbase*nbase/(nbase+1) + sd/(nbase+1);
            varbase = varbase*nbase/(nbase+1) + sd*sd/(nbase+1);
            sdbase = sqrt(varbase - avgbase*avgbase);
            nbase += 1;
        }
    }

    // remove points if offset precedes onset at the beginning
    if (ton[0] >= toff[0]){
        toff.remove(0);
    }
    if (ton.length() > toff.length()){
        ton.remove(ton.length()-1);
    }
    // delete sound whose duration is < minDur
    for(int ii=ton.length()-1; ii>=0; ii-=1){
        if (toff[ii] - ton[ii] < minDur/1000){
            ton.remove(ii);
            toff.remove(ii);
        }
    }

    // delete gap whose duration is < minISI
    for(int ii=ton.length()-1; ii>0; ii-=1){
        if (ton[ii] - toff[ii-1] < minISI/1000){
            ton.remove(ii);
            toff.remove(ii-1);
        }
    }

    QVector<QVector<double>> tmpout;
    tmpout.resize(2);
    tmpout[0].resize(ton.length());
    tmpout[1].resize(toff.length());
    for(int ii=0; ii<ton.length(); ii++){
        tmpout[0][ii] = ton[ii];
        tmpout[1][ii] = toff[ii];
    }

    widgetchart0->prepareSpecImage();

    //saveData2D(tmpout);
}

void MainWindow::openDir(QString pathname){
    QString filePathtmp;
    //QString filePathtmp2;

    if (!pathname.isEmpty()){
        filePathtmp = pathname;
    } else {
        QUrl fileurl = QFileDialog::getExistingDirectoryUrl(this, "Open dir", QUrl(wav.filePath), QFileDialog::ShowDirsOnly);
        filePathtmp = fileurl.toLocalFile();
        //filePathtmp = fileurl.path();
        //QMessageBox::information(this, "file path", filePathtmp);
        //QMessageBox::information(this, "file path", filePathtmp2);

        //filePathtmp.remove(0,1);
    }

    if (filePathtmp.isEmpty()){
        return;
    } else {
        listwidget->clear();
        appbusy = true;
        this->setCursor(Qt::WaitCursor);
        QDir filedir(filePathtmp);
        QStringList filters;
        filters << "*.wav";
        filedir.setNameFilters(filters);
        QStringList strlist = filedir.entryList(filters, QDir::NoFilter, QDir::Time);
        listwidget->addItems(strlist);
    }
    wav.filePath = filePathtmp;
    this->setWindowTitle("masaWavRead: " + wav.filePath);
    appbusy = false;
    this->setCursor(Qt::ArrowCursor);
}

void MainWindow::init(){
    ton.clear();
    toff.clear();
}


void MainWindow::calcSpectrogram(QVector<double> wavey, QVector<double> wavex, struct Spec& f){
    appbusy = true;
    this->setCursor(Qt::WaitCursor);

    f.mag.clear();
    f.phase.clear();
    f.y.clear();
    f.x.clear();
    f.mag_max = 0;
    f.dx = 0;
    f.dy = 0;

    int freq = v[f.ind].samplingfreq;

    //QMessageBox::information(this, "file path", QString::number(freq));
    //QMessageBox::information(this, "file path", QString::number(winwid));
    //QMessageBox::information(this, "file path", QString::number(padwid));
    //QMessageBox::information(this, "file path", QString::number(shiftwid));

    double stt = wavex[0];
    //QMessageBox::information(this, "!Done", QString::number(stt));
    double ent = wavex[wavex.length()-1];
    //QMessageBox::information(this, "!Done", QString::number(ent));
    int stpnt = 0; // starting point of calculation
    //QMessageBox::information(this, "!Done", QString::number(stpnt));
    int enpnt = int((ent-stt)*freq-1); // end point of calculation
    //QMessageBox::information(this, "!Done", QString::number(enpnt));
    if (stpnt < 0){
        stpnt = 0;
    }
    if (enpnt > wavex.length()-1){
        enpnt = wavex.length()-1;
    }

    int winpntorg = int(f.winwid*freq); // f.winwid*44100
    int padpnt = int(f.padwid*freq);
    int shiftpnt = int(f.shiftwid*freq); // this truncation will change the calc times

    double max_freq = 1*freq;
    double min_freq = 1/f.winwid;
    double padCoef = f.padwid/f.winwid;
    double dpnt = f.winwid/f.shiftwid;
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
    //QMessageBox::information(this, "!Done", QString::number(winpnt)); // 100

    int Ncalc = int((enpnt - stpnt - winpntorg)/shiftpnt)+1;
    //QMessageBox::information(this, "!Done", QString::number(Ncalc)); // 841
    //Spec_mag.resize(Ncalc);
    //Spec_phase.resize(Ncalc);
    //Spec_x.resize(Ncalc);

    f.mag.resize(Ncalc);
    f.phase.resize(Ncalc);
    f.x.resize(Ncalc);
    //QMessageBox::information(this, "Unable to save file", QString::number(f.mag.length())); // 3574
    //QMessageBox::information(this, "Unable to save file", QString::number(f.phase.length())); // 3574
    //QMessageBox::information(this, "Unable to save file", QString::number(f.x.length())); // 3574
    //f.mag.resize(Ncalc);
    //f.phase.resize(Ncalc);
    //f.x.resize(Ncalc);
    //QMessageBox::information(this, "Unable to save file", QString::number(f.mag.length())); // 3574
    //QMessageBox::information(this, "Unable to save file", QString::number(f.phase.length())); // 3574
    //QMessageBox::information(this, "Unable to save file", QString::number(f.x.length())); // 3574

    //QMessageBox::information(this, "Unable to save file", QString::number(int((enpnt - stpnt - winpnt)/shiftpnt))+1);

    //QMessageBox::information(this, "!Done", "ready");
    // Use FFTW library
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * uint(winpnt));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * uint(winpnt));

    // Use this if you do FFT for a fixed N data multiple times (it takes a few sec for init but then it uses the fastest algorithm)
    //fftw_plan p = fftw_plan_dft_1d(winpnt, in, out, FFTW_FORWARD, FFTW_MEASURE);
    // Or, use this if you do FFT for multiple data with different size (it doesn't need init time, but may not be optimal)
    fftw_plan p = fftw_plan_dft_1d(winpnt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    double hanning = 1;

    int sindex = 0;
    int stwin = stpnt;
    int enwin = stwin + winpnt;
    do {
        // Prepare data
        for(int ii=0; ii<winpnt; ii++){
            if (ii < (padpnt-winpntorg)/2 || (padpnt+winpntorg)/2 <= ii){
                in[ii][0] = 0; // real value, padding
            } else {
                //hanning = 0.5 * (1 - cos(2 * PI * (ii-(padpnt-winpntorg)/2) / (winpntorg+(padpnt-winpntorg)/2-1)));
                hanning = 0.5 * (1 - cos(2 * PI * (ii-(padpnt-winpntorg)/2) / (winpntorg-1)));
                // Hunning = 0.5 * (1 - cos(2 * pi * x / (numpnts(Hunning)-1)));
                in[ii][0] = wavey[stwin+ii-(padpnt-winpntorg)/2] * hanning; // real value
                //QMessageBox::information(this, "val", QString::number(stwin+ii-(padpnt-winpntorg)/2));
                //QMessageBox::information(this, "val", QString::number(in[ii][0]));
            }
            in[ii][1] = 0; // imaginary value
        }
        fftw_execute(p); // FFT (real output: out[ii][0]; imaginary output: out[ii][1])

        // Spec save
        //Spec_mag[sindex].resize(winpnt/2); // magnitude in [quantity peak]
        //Spec_phase[sindex].resize(winpnt/2); // phase in [radians] (-pi ~ pi)
        f.mag[sindex].resize(winpnt/2); // magnitude in [quantity peak]
        f.phase[sindex].resize(winpnt/2); // phase in [radians] (-pi ~ pi)
        //f.mag[sindex].resize(winpnt/2); // magnitude in [quantity peak]
        //f.phase[sindex].resize(winpnt/2); // phase in [radians] (-pi ~ pi)
        for(int ii=0; ii<winpnt/2; ii++){
            if (ii == 0 || ii == winpnt/2){ // or you can include the Nyquist (Nyquist is discarded in the current code)
                //Spec_mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt; // DC and Nyquist ([rms] magnitude is the same)
                f.mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt; // DC and Nyquist ([rms] magnitude is the same)
                //f.mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt; // DC and Nyquist ([rms] magnitude is the same)
                //FFTmaglog[ii] = 20 * log(10,FFTmag[ii])/; // powerspectrum
                //FFTpower[ii] = FFTmag[ii]*FFTmag[ii]; // powerspectrum
            } else {
                //Spec_mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt * 2; // pos+neg frequencies (further dividing with sqrt(2) will give you [rms] magnitude)
                f.mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt * 2; // pos+neg frequencies (further dividing with sqrt(2) will give you [rms] magnitude)
                //f.mag[sindex][ii] = sqrt(out[ii][0]*out[ii][0] + out[ii][1]*out[ii][1])/winpnt * 2; // pos+neg frequencies (further dividing with sqrt(2) will give you [rms] magnitude)
                //FFTpower[ii] = FFTmag[ii]*FFTmag[ii]/2; // powerspectrum
            }
            //if (Spec_mag_max < Spec_mag[sindex][ii]){
            if (f.mag_max < f.mag[sindex][ii]){
                //Spec_mag_max = Spec_mag[sindex][ii];
                f.mag_max = f.mag[sindex][ii];
                //f.mag_max = f.mag[sindex][ii];
            }
            //Spec_phase[sindex][ii] = atan(out[ii][1] / out[ii][0]); // DC and Nyquist
            f.phase[sindex][ii] = atan(out[ii][1] / out[ii][0]); // DC and Nyquist
            //f.phase[sindex][ii] = atan(out[ii][1] / out[ii][0]); // DC and Nyquist
            //double bp = out[ii][1] / out[ii][0];
            //QMessageBox::information(this, "!Done", QString::number(bp));
            //bp = atan(bp);
            //QMessageBox::information(this, "!Done", QString::number(bp));
        }
        //Spec_x[sindex] = stt + f.winwid/2 + double(sindex * shiftpnt) / freq;
        f.x[sindex] = stt + f.winwid/2 + double(sindex * shiftpnt) / freq;
        //f.x[sindex] = stt + f.winwid/2 + double(sindex * shiftpnt) / freq;

        sindex++;
        stwin = stpnt + sindex * shiftpnt;
        enwin = stwin + winpnt;
    } while(sindex < Ncalc); // Ncalc rather than Ncalc-1

    //Spec_y.resize(winpnt/2);
    f.y.resize(winpnt/2);
    //f.y.resize(winpnt/2);
    for(int ii=0; ii<winpnt/2; ii++){
        //Spec_y[ii] = 0 + ii*df;
        f.y[ii] = 0 + ii*df;
        //f.y[ii] = 0 + ii*df;
    }
    //Spec_dx = f.shiftwid;
    f.dx = f.shiftwid;
    //QMessageBox::information(this, "!Done", QString::number(Spec_dx));
    //Spec_dy = df;
    f.dy = df;
    //QMessageBox::information(this, "!Done",  QString::number(Spec_dy));
    //QMessageBox::information(this, "!Done", QString::number(Spec_x.length()));
    //QMessageBox::information(this, "!Done", QString::number(Spec_y.length()));
    //f.dx = f.shiftwid;
    //f.dy = df;

    // Destroy FFTW
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    appbusy = false;
    this->setCursor(Qt::ArrowCursor);

}


void MainWindow::saveData1D(QVector<double> vec, QString fileName){
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    //QString rhdfilename = parentwindow->getFileName();
    QString rhdfilename = getFileName();
    QString rhdpathname = getFilePathName();
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

        for(int ii=0; ii<vec.length(); ii++){
            out << vec[ii] << "\r\n";
        }

        file.close();
    }
}

void MainWindow::saveData2D(QVector<QVector<double>> vec, QString fileName){
    //parentwindow = (MainWindow*)parentWidget(); // this makes it crash, but I don't know why
    //QString rhdfilename = parentwindow->getFileName();
    QString rhdfilename = getFileName();
    QString rhdpathname = getFilePathName();
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

        int maxlen = 0;
        for(int ii=0; ii<vec.length(); ii++){
            if (maxlen < vec[ii].length()){
                maxlen = vec[ii].length();
            }
        }
        for(int jj=0; jj<maxlen; jj++){
            for(int ii=0; ii<vec.length(); ii++){
                out << vec[ii][jj]; // real value
                if (jj<vec[ii].length()-1){
                    out << "\t";
                }
            }
            out << "\r\n";
        }

        file.close();
    }
}


void MainWindow::pushPlaySound(){
    if (appbusy || wav.filePath.isEmpty()){
        return;
    }

    QString fileNametmp = listwidget->currentItem()->text();
    QString filePathtmp = wav.filePath + "/" + fileNametmp;

    player.setAudioRole(QAudio::MusicRole);
    player.setMedia(QUrl::fromLocalFile(filePathtmp));
    player.setPlaybackRate(1.0);
    player.setVolume(100);
    QObject::connect(&player, SIGNAL(mediaStatusChanged(QMediaPlayer::MediaStatus)), this, SLOT(playClicked()));

    /* Also see
     * https://www.ipentec.com/document/csharp-play-wave-file-using-wave-api
     * https://docs.microsoft.com/en-us/windows/desktop/directshow/audio-streaming-sample-code
    // https://doc.qt.io/qt-5.11/qtmultimedia-multimediawidgets-player-example.html
     */
}

void MainWindow::calcSnipWAV(){
    double stt = widgetchart0->getXmin();
    double ent = widgetchart0->getXmax();
    int stpnt = int(stt * wav.SamplesPerSec);
    int enpnt = int(ent * wav.SamplesPerSec);
    if (stpnt > enpnt){
        QMessageBox::information(this, "!Error", "Out of range");
        return;
    }
    if (stpnt < 0){
        stpnt = 0;
    }
    if (enpnt >= wav.dataL.length()){
        enpnt = wav.dataL.length() - 1;
    }
    int lenpnt = enpnt - stpnt + 1;  // including enpnt
    wavSnipL.resize(lenpnt);
    if (wav.NCh == 1){
        for(int ii=0; ii<lenpnt; ii++){
            wavSnipL[ii] = wav.dataL[ii+stpnt]; // 2
        }
    } else if (wav.NCh == 2){
        wavSnipR.resize(lenpnt);
        for(int ii=0; ii<lenpnt; ii++){
            wavSnipL[ii] = wav.dataL[ii+stpnt]; // 2
            wavSnipR[ii] = wav.dataR[ii+stpnt]; // 2
        }
    }
}

void MainWindow::pushPlaySoundPart(){
    if (appbusy || wav.filePath.isEmpty()){
        return;
    }
    player.pause();
    player.setAudioRole(QAudio::MusicRole);

    calcSnipWAV();
    int lenpnt = wavSnipL.length();
    if (lenpnt <= 0) return;
    qint32 chunksize = lenpnt * (wav.NCh * wav.bitsPerSample/8) - 46; // *** This determines the duration of playback. I don't know why but -46 (header size) works the best so far, although this is less than the actual data size.
    qint32 wavsize = lenpnt + 36; // SAP data format. Not sure if this is general wav format. I think this does not affect the playback anyway.

    buffer = new QBuffer();
    buffer->open(QIODevice::ReadWrite);
    buffer->seek(0);

    // write parametes
    QDataStream out(buffer);
    out.setVersion(QDataStream::Qt_4_8); // should I use pointer? in -> setVersion ?
    out.setByteOrder(QDataStream::LittleEndian);
    out.setFloatingPointPrecision(QDataStream::SinglePrecision);
    out << wav.RIFF[0] << wav.RIFF[1] << wav.RIFF[2] << wav.RIFF[3]; //4
    out << wavsize; // 4
    out << wav.WAVE[0] << wav.WAVE[1] << wav.WAVE[2] << wav.WAVE[3]; // 4
    out << wav.fmt[0] << wav.fmt[1] << wav.fmt[2] << wav.fmt[3]; // 4
    out << wav.ChunkSize; // 4
    out << wav.AudioFormat ; // 2
    out << wav.NCh; // 2
    out << wav.SamplesPerSec; // 4
    out << wav.bytesPerSec; // 4
    out << wav.blockAlign ; // 2
    out << wav.bitsPerSample; // 2
    out << wav.Subchunk2ID[0] << wav.Subchunk2ID[1] << wav.Subchunk2ID[2] << wav.Subchunk2ID[3]; // 4
    out << chunksize; // 4
    if (wav.NCh == 1){
        for(int ii=0; ii<lenpnt; ii++){
            out << wavSnipL[ii]; // 2
        }
    } else if (wav.NCh == 2){
        for(int ii=0; ii<lenpnt; ii++){
            out << wavSnipL[ii];
            out << wavSnipR[ii];
        }
    }

    buffer->seek(0);
    //QMessageBox::information(this, "!Done", QString::number(stt) + "s[" + QString::number(stpnt) + "] - " + QString::number(ent) + "s[" + QString::number(enpnt) + "]" + "\r\nbuffer size: " + QString::number(buffer->size()));

    player.setPlaybackRate(1.0);
    player.setVolume(100);
    player.setMedia(QMediaContent(), buffer); // buffer needs a wav header
    QObject::connect(&player, SIGNAL(mediaStatusChanged(QMediaPlayer::MediaStatus)), this, SLOT(playClicked()));
    //QMessageBox::information(this, "!Done", QString::number(player.mediaStatus()));
}

void MainWindow::playClicked(){
    if (player.mediaStatus() == QMediaPlayer::LoadedMedia){
        player.play();
    }
}

void MainWindow::pushStopSound(){
    player.pause(); // if stopped, status changes and thus playClicked() is called to play sound again
}

void MainWindow::pushSaveWAV(){

    QString rhdfilename = getFileName();
    QString rhdpathname = getFilePathName();
    rhdfilename.chop(4); // remove ".rhd"
    QString tmpfile;

    calcSnipWAV();
    int lenpnt = wavSnipL.length();
    if (lenpnt <= 0) return;
    qint32 chunksize = lenpnt * (wav.NCh * wav.bitsPerSample/8); // *** Although I add -46 for playnback, I don't do that here as this is the actual data size.
    qint32 wavsize = lenpnt + 36; // SAP data format. Not sure if this is general wav format. It makes more sense if I add 38?

    buffer = new QBuffer();
    buffer->open(QIODevice::ReadWrite);
    buffer->seek(0);

    if (tmpfile.isEmpty()){
        tmpfile = QFileDialog::getSaveFileName(this, "Save as:", rhdpathname + "/" + rhdfilename + "_test.wav", "ASCII (*.wav);;All Files (*)");
    }
    rhdfilename= rhdfilename + ".wav";
    if (tmpfile.isEmpty()){
        return;
    } else {
        QFile file(tmpfile);
        if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, "Unable to save file", file.errorString());
            return;
        }
        // write parametes
        QDataStream out(&file);
        out.setVersion(QDataStream::Qt_4_8); // should I use pointer? in -> setVersion ?
        out.setByteOrder(QDataStream::LittleEndian);
        out.setFloatingPointPrecision(QDataStream::SinglePrecision);
        out << wav.RIFF[0] << wav.RIFF[1] << wav.RIFF[2] << wav.RIFF[3]; //4
        out << wavsize; // 4
        out << wav.WAVE[0] << wav.WAVE[1] << wav.WAVE[2] << wav.WAVE[3]; // 4
        out << wav.fmt[0] << wav.fmt[1] << wav.fmt[2] << wav.fmt[3]; // 4
        out << wav.ChunkSize; // 4
        out << wav.AudioFormat ; // 2
        out << wav.NCh; // 2
        out << wav.SamplesPerSec; // 4
        out << wav.bytesPerSec; // 4
        out << wav.blockAlign ; // 2
        out << wav.bitsPerSample; // 2
        out << wav.Subchunk2ID[0] << wav.Subchunk2ID[1] << wav.Subchunk2ID[2] << wav.Subchunk2ID[3]; // 4
        out << chunksize; // 4
        if (wav.NCh == 1){
            for(int ii=0; ii<wavSnipL.size(); ii++){
                out << wavSnipL[ii]; // 2
            }
        } else if (wav.NCh == 2){
            for(int ii=0; ii<wavSnipL.size(); ii++){
                out << wavSnipL[ii];
                out << wavSnipR[ii];
            }
        }
        file.close();
    }
}

QVector<qint16>& MainWindow::getWavraw(){
    return wav.dataL;
}

QVector<double>& MainWindow::getDatat(){
    return wav.wavt;
}

QVector<double>& MainWindow::getSeries(int a){
    return v[a].data;
}

QVector<double>& MainWindow::getONt(){
    return ton;
}
QVector<double>& MainWindow::getOFFt(){
    return toff;
}

bool MainWindow::getFileOpened(){
    return fileopened;
}

QString MainWindow::getFileName(){
    return wav.fileName;
}

QString MainWindow::getFilePathName(){
    return wav.filePath;
}
struct Spec& MainWindow::getSpec(){
//    return Spec_x;
    return specs[0];
}
QVector<double>& MainWindow::getSpec_x(){
    //return Spec_x;
    return specs[0].x;
}
QVector<double>& MainWindow::getSpec_y(){
    //return Spec_y;
    return specs[0].y;
}
QVector<QVector<double>>& MainWindow::getSpec_mag(){
    //return Spec_mag;
    return specs[0].mag;
}
QVector<QVector<double>>& MainWindow::getSpec_phase(){
    //return Spec_phase;
    return specs[0].phase;
}
double MainWindow::getSpec_dx(){
    //return Spec_dx;
    return specs[0].dx;
}
double MainWindow::getSpec_mag_max(){
    //return Spec_mag_max;
    return specs[0].mag_max;
}
double MainWindow::getSpec_dy(){
    //return Spec_dy;
    return specs[0].dy;
}

QVector<QString> MainWindow::myDialog(QString dialogtitle, QVector<QString> dialoglabel, QVector<QString> dialogvalue){
    QVector<QString> returnvalue;
    if (dialoglabel.length() != dialogvalue.length()){
        return returnvalue;
    }
    int nfield = dialoglabel.length();

    QDialog dialog(this);
    dialog.setWindowTitle(dialogtitle);
    QFormLayout form(&dialog);
    //form.addRow(new QLabel("Parameters:"));
    QList<QLineEdit*> fields;
    QLineEdit* lineEdit[nfield];
    for(int ii = 0; ii < nfield; ++ii) {
        lineEdit[ii] = new QLineEdit(dialogvalue[ii], &dialog);
        form.addRow(dialoglabel[ii], lineEdit[ii]);
        fields << lineEdit[ii];
    }
    QDialogButtonBox buttonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dialog);
    form.addRow(&buttonBox);
    QObject::connect(&buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    QObject::connect(&buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

    if (dialog.exec() != QDialog::Accepted) {
        return returnvalue;
    } else {
        returnvalue.resize(nfield);
        for(int ii = 0; ii < nfield; ++ii) {
            returnvalue[ii] = lineEdit[ii]->text();
        }
        return returnvalue;
    }

}

void MainWindow::specsliderChanged(){
    widgetchart0->reprint();
}

int MainWindow::getSpecslider(){
    return specslider->value();
}

