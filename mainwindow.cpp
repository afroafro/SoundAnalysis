#include "QWidgetWAVChart.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
//#include "QChartViewWAV.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    // parameters (these may not be working; look at WidgetIntanChart.cpp)
    fileopened = false;
    appbusy = false;
    /*/////////////////////////////////////*/

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

    QWidget* widgetchartbuttons = new QWidget(this);

    // chart buttons
    QGridLayout *chartButtonLayout = new QGridLayout;
    const int nChartButton = 8;
    QPushButton *chartButtons[nChartButton];
    chartButtons[0] = new QPushButton("top");
    QObject::connect(chartButtons[0], SIGNAL(released()), widgetchart0, SLOT(releasedChartButtonTop()));
    chartButtons[1] = new QPushButton("up");
    QObject::connect(chartButtons[1], SIGNAL(released()), widgetchart0, SLOT(releasedChartButtonUp()));
    chartButtons[2] = new QPushButton("down");
    QObject::connect(chartButtons[2], SIGNAL(released()), widgetchart0, SLOT(releasedChartButtonDown()));
    chartButtons[3] = new QPushButton("bottom");
    QObject::connect(chartButtons[3], SIGNAL(released()), widgetchart0, SLOT(releasedChartButtonBottom()));
    chartButtons[4] = new QPushButton("4");
    QObject::connect(chartButtons[4], SIGNAL(released()), widgetchart0, SLOT(releasedChartButton4()));
    chartButtons[5] = new QPushButton("8");
    QObject::connect(chartButtons[5], SIGNAL(released()), widgetchart0, SLOT(releasedChartButton8()));
    chartButtons[6] = new QPushButton("12");
    QObject::connect(chartButtons[6], SIGNAL(released()), widgetchart0, SLOT(releasedChartButton12()));
    chartButtons[7] = new QPushButton("16");
    QObject::connect(chartButtons[7], SIGNAL(released()), widgetchart0, SLOT(releasedChartButton16()));
    for (int ii=0; ii<nChartButton; ++ii){
        chartButtonLayout->addWidget(chartButtons[ii], ii, 0);
    }

    widgetchartbuttons->setLayout(chartButtonLayout);

    mainLayout->addWidget(widgetchartbuttons, 0, 0);
    mainLayout->addWidget(widgetchart0, 0, 1);

    // *** initialize widget
    //widgetchart0->initParams();
}

MainWindow::~MainWindow() {
    delete ui;
}


void MainWindow::pushFile(QAction *action){
    int value = action->data().toInt();
    if (value == 0){
        openFile();
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

void MainWindow::pushDisplay(QAction *action){
    int value = action->data().toInt();
    if (value == 0){
    } else if (value == 1){
    } else if (value == 2){ // show parameters
    } else if (value == 3){ // show sorted spikes
    }
}


void MainWindow::openFile(QString rhdfilename, QString rhdpathname){
    //QMessageBox::information(this, "!!", "pressed");
    if (appbusy){
        return;
    }

    QString fileNametmp;
    QString filePathtmp;

    if (!rhdfilename.isEmpty() && !rhdpathname.isEmpty()){
        fileNametmp = rhdfilename;
        filePathtmp = rhdpathname;
    } else {
        QUrl fileurl = QFileDialog::getOpenFileUrl(this, "Open", QUrl(filePath), "wav file (*.wav);;All Files (*)");
        fileNametmp = fileurl.fileName();
        filePathtmp = fileurl.path();
        filePathtmp.remove(0,1);
    }
    //QMessageBox::information(this, "Unable to open file", fileName);
    //QMessageBox::information(this, "Unable to open file", filePath);

    if (fileNametmp.isEmpty()){
        return;
    } else {
        appbusy = true;
        this->setCursor(Qt::WaitCursor);
        QFile file(filePathtmp);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(this, "Unable to open file", file.errorString());
            return;
        }

        QDataStream in(&file);
        in.setVersion(QDataStream::Qt_4_8); // should I use pointer? in -> setVersion ?
        in.setByteOrder(QDataStream::LittleEndian);
        in.setFloatingPointPrecision(QDataStream::SinglePrecision);

        filebytesize = file.size();

        in >> RIFF[0] >> RIFF[1] >> RIFF[2] >> RIFF[3];
        in >> wavfilesize;
        in >> WAVE[0] >> WAVE[1] >> WAVE[2] >> WAVE[3];
        in >> fmt[0] >> fmt[1] >> fmt[2] >> fmt[3];
        in >> ChunkSize; // 16
        in >> AudioFormat ;
        in >> NCh;
        in >> SamplesPerSec;
        in >> bytesPerSec;
        in >> blockAlign ;
        in >> bitsPerSample;
        in >> Subchunk2ID[0] >> Subchunk2ID[1] >> Subchunk2ID[2] >> Subchunk2ID[3];
        in >> Subchunk2Size;
        if (NCh == 1){
            wavDataL.resize(Subchunk2Size/2);
            for(int ii=0; ii<wavDataL.size(); ii++){
                in >> wavDataL[ii];
            }
        } else if (NCh == 2){
            wavDataL.resize(Subchunk2Size/4);
            wavDataR.resize(Subchunk2Size/4);
            for(int ii=0; ii<wavDataL.size(); ii++){
                in >> wavDataL[ii];
                in >> wavDataR[ii];
            }
        } else {
            // error
            file.close();
            QMessageBox::information(this, "!Exception", "Seems not match any data format.");
            appbusy = false;
            return;
        }
        file.close();

        QString msg = "";
        msg += "RIFF:\t" + QString::number(RIFF[0]) + QString::number(RIFF[1]) + QString::number(RIFF[2]) + QString::number(RIFF[3]) + "\n";
        msg += "wavfile:\t" + QString::number(wavfilesize) + "\n";
        msg += "WAVE:\t" + QString::number(WAVE[0]) + QString::number(WAVE[1]) + QString::number(WAVE[2]) + QString::number(WAVE[3]) + "\n";
        msg += "fmt:\t" + QString::number(fmt[0]) + QString::number(fmt[1]) + QString::number(fmt[2]) + QString::number(fmt[3]) + "\n";
        msg += "Chunksize:\t" + QString::number(ChunkSize) + "\n";
        msg += "audio:\t" + QString::number(AudioFormat) + "\n";
        msg += "NCh:\t" + QString::number(NCh) + "\n";
        msg += "samples/s:\t" + QString::number(SamplesPerSec) + "\n";
        msg += "bytes/s:\t" + QString::number(bytesPerSec) + "\n";
        msg += "blockalign:\t" + QString::number(blockAlign) + "\n";
        msg += "bits/sample:\t" + QString::number(bitsPerSample) + "\n";
        msg += "subchunk:\t" + QString::number(Subchunk2ID[0]) + QString::number(Subchunk2ID[1]) + QString::number(Subchunk2ID[2]) + QString::number(Subchunk2ID[3]) + "\n";
        msg += "subchunksize:\t" + QString::number(Subchunk2Size) + "\n";
        QMessageBox::information(this, "!Done", msg);
    }

    fileName = fileNametmp;
    filePath = filePathtmp;
    QStringList dirnames = filePathtmp.split("/");
    fileDir = "";
    for (int ii=0; ii<dirnames.length()-1; ++ii){
        fileDir = fileDir + dirnames[ii] + "/";
    }

    fileopened = true;

    wavt.resize(wavDataL.length());
    for (int ii=0; ii<wavDataL.length(); ii++){
        wavt[ii] = double(ii)/double(SamplesPerSec);
    }

    this->setWindowTitle("masaWavRead: " + fileName);
    appbusy = false;
    this->setCursor(Qt::ArrowCursor);

    widgetchart0->setData();
}

QVector<qint16>& MainWindow::getWavraw(){
    return wavDataL;
}

QVector<double>& MainWindow::getDatat(){
    return wavt;
}

bool MainWindow::getFileOpened(){
    return fileopened;
}

QString MainWindow::getFileName(){
    return fileName;
}

QString MainWindow::getFilePathName(){
    return filePath;
}
