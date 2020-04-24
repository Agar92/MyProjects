#include <QVBoxLayout>
#include <QLineEdit>
#include <QPushButton>

#include "AWidget.h"

AWidget::AWidget(QWidget *parent) :
        QWidget(parent)
{

    QVBoxLayout *layout = new QVBoxLayout(this);
    //layout->setSpacing(0);

    time_widget = new QLineEdit(this);
    time_widget->setFont(QFont("Times", 22, QFont::Bold));
    time_widget->setStyleSheet("color: blue;");
    time_widget->setAlignment(Qt::AlignCenter);


    timer = new QTimer();
    QObject::connect(timer, SIGNAL(timeout()), this, SLOT(time()));

    QPushButton *Start = new QPushButton("Start", this);
    Start->setFont(QFont("Times", 22, QFont::Bold));
    //Start->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    QPushButton *Pause = new QPushButton("Pause", this);
    Pause->setFont(QFont("Times", 22, QFont::Bold));
    //Pause->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    QPushButton *Stop = new QPushButton("Stop", this);
    Stop->setFont(QFont("Times", 22, QFont::Bold));
    QPushButton *Restart = new QPushButton("Restart", this);
    Restart->setFont(QFont("Times", 22, QFont::Bold));
    Restart->setStyleSheet("color: red;");

    QObject::connect(Start, SIGNAL(clicked()), this, SLOT(start()));
    QObject::connect(Pause, SIGNAL(clicked()), this, SLOT(pause()));
    QObject::connect(Stop, SIGNAL(clicked()), this, SLOT(stop()));
    QObject::connect(Restart, SIGNAL(clicked()), this, SLOT(restart()));

    layout->addWidget(time_widget);
    layout->addWidget(Start);
    layout->addWidget(Pause);
    layout->addWidget(Stop);
    layout->addWidget(Restart);

    this->setStyleSheet("color: blue");
    this->setLayout(layout);
    this->setWindowTitle("Секундомер");
    this->resize(200, 100);
}

void AWidget::start(){
   if(pause_stop==1){
       timer->start(1000);
    }
    if(pause_stop==2){
        time_widget->setStyleSheet("color: red;");
    }
    if(pause_stop==0){
     timer->start(1000);
    }
}

void AWidget::pause()
{
    timer->stop();
    pause_stop=1;
}

void AWidget::stop()
{
    timer->stop();
    pause_stop=2;
}

void AWidget::restart()
{
    pause_stop=0;
}

void AWidget::time()
{
    if(++s == 60)
    {
        s = 0;
        if(++m == 60)
        {
            m = 0;
            ++h;
        }
    }
    time_widget->setText(QString("%1:%2:%3")
        .arg(h, 2, 10, QLatin1Char('0'))
        .arg(m, 2, 10, QLatin1Char('0'))
        .arg(s, 2, 10, QLatin1Char('0')));
}
