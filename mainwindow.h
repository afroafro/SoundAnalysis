#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <QMediaPlayer>

class QWidgetWAVChart; // forward declaration
class windows;
class IMultiMediaStream;

struct WAVdata { // data read from .wav files
    QString filePath;
    QString fileName;
    QString fileDir;
    qint64 filebytesize;
    qint8    RIFF[4];        // RIFF Header Magic header
    qint32   wavfilesize;    // 4 + (8 + SubChunk1Size) + (8 + SubChunk2Size), thus 36 + SubChunk2Size, according to http://soundfile.sapp.org/doc/WaveFormat/. However, in SAP data it is numpnts(data) + 36
    qint8    WAVE[4];        // WAVE Header
    /* "fmt" sub-chunk */
    qint8    fmt[4];         // FMT header
    qint32   ChunkSize;      // RIFF Chunk Size (length of above = 16 for PCM)
    qint16   AudioFormat;    // Audio format 1=PCM (i.e. Linear quantization), other numbers: compression (6=mulaw,7=alaw,257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM)
    qint16   NCh;      // Number of channels 1=Mono 2=Sterio
    qint32   SamplesPerSec;  // Sampling Frequency in Hz
    qint32   bytesPerSec;    // bytes per second == SampleRate * NumChannels * BitsPerSample/8
    qint16   blockAlign;     // 2=16-bit mono, 4=16-bit stereo == NumChannels * BitsPerSample/8
    qint16   bitsPerSample;  // Number of bits per sample 8 bits = 8, 16 bits = 16, etc.
    /* "data" sub-chunk */
    qint8    Subchunk2ID[4]; // "data"  string
    qint32   Subchunk2Size;  // Sampled data length == NumSamples * NumChannels * BitsPerSample/8
    int datasamplesize = 0;
    QVector<qint16> dataL;
    QVector<qint16> dataR;
    QVector<double> wavt; // Do not calc from dx*xpnt (as dx is rounded), but calc xpnt/samplerate
};
struct MyVector { // parameters and results of song analysis
    QString name;
    QVector<double> data;
    int samplingfreq;
    int xind = 0; // this defines x wave
    int yaxisind = -1; // Target axis (if -1, exclude from graphs. if >= 0, show in AxisY[yaxisind])
};
struct MyXVector { // parameters and results of song analysis
    QString name;
    QVector<double> data;
};
struct SongPart { //
    QVector<double> ton;
    QVector<double> toff;
};
struct Spec { // parameters and results of spectrogram
    // parameters
    int ind = 0; // target index (mydata[ind]), freq = v[f.ind].samplingfreq;
    double winwid = 0.01;
    double padwid = 0.01;
    double shiftwid = 0.005;
    // data
    QVector<double> x;
    QVector<double> y;
    QVector<QVector<double>> mag;
    QVector<QVector<double>> phase;
    double dx = 0;
    double dy = 0;
    double mag_max = 0;
};
struct FiltParam { // parameters for filtering
    double samplerate = 44100;
    double HPfreq = 200;
    double HPattn = 50;
    double HPripple = 0.1;
    double LPfreq = 20000;
    double LPattn = 50;
    double LPripple = 0.1;
};

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
    public:
        explicit MainWindow(QWidget *parent = nullptr);
        ~MainWindow();
        void openDir(QString="");
        //void openFile(QString="");
        bool getFileOpened();
        QVector<qint16>& getWavraw();
        QVector<double>& getDatat();
        QVector<double>& getONt();
        QVector<double>& getOFFt();
        QVector<double>& getSeries(int);
        QVector<double>& getSpec_x();
        QVector<double>& getSpec_y();
        QVector<QVector<double>>& getSpec_mag();
        QVector<QVector<double>>& getSpec_phase();
        double getSpec_dx();
        double getSpec_dy();
        struct Spec& getSpec();
        double getSpec_mag_max();
        bool appbusy;

        QString getFileName();
        QString getFilePathName();
        void saveData1D(QVector<double>, QString="");
        void saveData2D(QVector<QVector<double>>, QString="");
        int getSpecslider();

private:
       Ui::MainWindow *ui;
       QWidget* widget;
       QWidget* buttonwidget;
       QGridLayout* mainLayout;
       QGridLayout* bottomLayout;
       QWidgetWAVChart* widgetchart0;
       QListWidget* listwidget;
       QWidget* bottomwidget;
       QVector<QSignalMapper*> buttonsignalmapper; // used to be [16][3] array, but changed to a vector
       QSlider* specslider;
       QVector<QPushButton*> buttons;

       bool fileopened;
       bool userstopped;

       struct WAVdata wav;
       QVector<struct Spec> specs;
       QVector<struct MyVector> v;
       QVector<struct MyXVector> x;

       QVector<QVector<double>> anaDataL;
       QVector<QVector<double>> anaDataR;
       QVector<qint16> wavSnipL; // for calculation of the shown figure
       QVector<qint16> wavSnipR; // for calculation of the shown figure
       QVector<double> ton;
       QVector<double> toff;

       const double PI = 3.141592653589793238463;

       // for default parameters (necessary?)
       double highpass_freq = 200;
       double highpass_attn = 50;
       double highpass_ripple = 0.1;
       double lowpass_freq = 15000;
       double lowpass_attn = 50;
       double lowpass_ripple = 0.1;

       void enableButtons();
       void disableButtons();
       void calcSnipWAV();
       void init();
       void runBatch();
       void checkBatch();
       void addAnaWave(QVector<double>);
       QVector<double> filterWave(QVector<double>, struct FiltParam);
       QVector<double> absWave(QVector<double>);
       QVector<double> differentiateWave(QVector<double>);
       QVector<double> avgDiscretizeWave(QVector<double>, int);
       QVector<QString> myDialog(QString, QVector<QString>, QVector<QString>);
       void calcSpectrogram(QVector<double>, QVector<double>, struct Spec&);

       QMediaPlayer player;
       QBuffer* buffer;

    public slots: // public, private are ignored for Qt slots. All slots are actually public
        void pushFile(QAction*);
        void pushDisplay(QAction*);
        void pushAmp();
        void pushFilter();
        void pushWaves();
        void openFile();
        void pushBatch();
        void pushStopBatch();
        void detectSound();
        void pushPlaySound();
        void pushStopSound();
        void playClicked();
        void pushPlaySoundPart();
        void pushSaveWAV();
        void pushRhythm();
        void specsliderChanged();
};

#endif // MAINWINDOW_H
