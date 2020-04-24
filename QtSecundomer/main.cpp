#include "mainwindow.h"
#include <AWidget.h>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.resize(300,200);
    w.show();
    return a.exec();
}
