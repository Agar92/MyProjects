#ifndef AWIDGET_H
#define AWIDGET_H

#include <QLineEdit>
#include <QTime>
#include <QTimer>

class AWidget : public QWidget
{
    Q_OBJECT

public:
    explicit AWidget(QWidget *parent = 0);

public slots:
    void start();
    void restart();
    void stop();
    void pause();
    void time();

private:
    QLineEdit *time_widget;
    QTimer *timer;
    int pause_stop=0;
    int h;
    int m;
    int s;
};

#endif // AWIDGET_H
