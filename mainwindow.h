#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>

class QWidgetWAVChart; // forward declaration

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
    public:
        explicit MainWindow(QWidget *parent = nullptr);
        ~MainWindow();
        void openFile(QString="", QString="");
        bool getFileOpened();
        QVector<qint16>& getWavraw();
        QVector<double>& getDatat();
        QString getFileName();
        QString getFilePathName();

private:
       Ui::MainWindow *ui;
       QWidget* widget;
       QGridLayout* mainLayout;
       QWidgetWAVChart* widgetchart0;

       QString filePath;
       QString fileName;
       QString fileDir;
       qint64 filebytesize;
       bool fileopened;
       bool appbusy;

       qint8    RIFF[4];        // RIFF Header Magic header
       qint32   wavfilesize;
       qint8    WAVE[4];        // WAVE Header
       /* "fmt" sub-chunk */
       qint8    fmt[4];         // FMT header
       qint32   ChunkSize;      // RIFF Chunk Size (length of above = 16)
       qint16   AudioFormat;    // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
       qint16   NCh;      // Number of channels 1=Mono 2=Sterio
       qint32   SamplesPerSec;  // Sampling Frequency in Hz
       qint32   bytesPerSec;    // bytes per second
       qint16   blockAlign;     // 2=16-bit mono, 4=16-bit stereo
       qint16   bitsPerSample;  // Number of bits per sample
       /* "data" sub-chunk */
       qint8    Subchunk2ID[4]; // "data"  string
       qint32   Subchunk2Size;  // Sampled data length
       QVector<qint16> wavDataL;
       QVector<qint16> wavDataR;
       QVector<double> wavt;

    public slots: // public, private are ignored for Qt slots. All slots are actually public
        void pushFile(QAction*);
        void pushDisplay(QAction*);
};

#endif // MAINWINDOW_H
