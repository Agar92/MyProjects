#include "mainwindow.h"

#include <QApplication>

//https://qtcentre.org/threads/60272-Signals-and-Slots-across-different-classes-(Qt5-Version)
//Example of complex connect() usage.

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
