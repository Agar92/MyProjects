#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMap>
#include <QList>
#include <QMediaPlayer>
#include <QMenu>
#include <QAction>
#include <QTimer>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    void keyPressEvent(QKeyEvent *e);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_horizontalSlider_actionTriggered(int action);

    void on_quitButton_clicked();

    void on_actionNew_word_triggered();

    void on_actionExit_triggered();

    void on_timeEdit_userTimeChanged(const QTime &time);

    void slotTimeout();

private:
    QAction * newAct;
    Ui::MainWindow *ui;
    QMap<QString,QPixmap> Map;
    QList<QString> names;
    int curnum;
    int size;
    QTimer * timer;
    QMediaPlayer * music;
};

#endif // MAINWINDOW_H
