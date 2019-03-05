#include "./qlistwidgetwav.h"

QListWidgetWAV::QListWidgetWAV(QWidget *parent) : QListWidget(parent){
//    parentwindow = (MainWindow*)parentWidget();
//    QObject::connect(this, SIGNAL(itemSelectionChanged()), this, SLOT(callOpenFile()));
}
/*
void QListWidgetWAV::itemSelectionChanged(){
    QListWidget::itemSelectionChanged();
}

void QListWidgetWAV::callOpenFile(){
    QMessageBox::information(this, "!!", "pressed");
    QList<QListWidgetItem*> list;
    list = this->selectedItems();
    if (list.length()>0){
        parentwindow->openFile(list[0]->text());
    }
}
*/
