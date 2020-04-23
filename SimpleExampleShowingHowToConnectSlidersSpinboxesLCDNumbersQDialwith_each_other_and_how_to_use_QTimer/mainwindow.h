#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QTime>
#include <QTimer>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void ConnectSignals();
    void DisconnectSignals();

private slots:
    void on_ViewButton_clicked();

private slots:
    void Valuechanged(int val);
    void UpdateTime();

    void on_dial_valueChanged(int value);

private:
    Ui::MainWindow *ui;
    QTimer * timer;
};

#endif // MAINWINDOW_H
