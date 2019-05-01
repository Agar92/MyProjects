#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPixmap>
#include <QMap>
#include <QList>
#include <QIcon>
#include <stdlib.h>
#include <QDebug>
#include <QMessageBox>
#include <QKeyEvent>
#include <QMediaPlayer>
#include <QFileDialog>
#include <QUrl>
#include <QInputDialog>
#include <QTimer>
#include <QTime>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    //our pics:
    QPixmap pix_hare(":/images/hare.jpg");
    QPixmap pix_sun(":/images/sun.jpg");
    QPixmap pix_train(":/images/train.png");
    QPixmap pix_tree(":/images/tree.png");
    //insert:
    Map.insert("hare",pix_hare);
    Map.insert("sun",pix_sun);
    Map.insert("train",pix_train);
    Map.insert("tree",pix_tree);
    //names:
    names.push_back("hare");
    names.push_back("sun");
    names.push_back("train");
    names.push_back("tree");
    size=Map.size();
    ui->setupUi(this);
    curnum=rand()%size;
    QString name=names[curnum];
    //qDebug()<<"curnum="<<curnum;
    //curnum=2;
    ui->label->setPixmap(Map[name]);
    QIcon quit(":/images/quit.png");
    ui->quitButton->setIcon(quit);
    music = new QMediaPlayer();
    music->setMedia(QUrl("qrc:/sounds/river_flows_in_you.mp3"));
    music->setVolume(20);
    ui->volumeSlider->setValue(20);
    connect(ui->volumeSlider, SIGNAL(valueChanged(int)),music, SLOT(setVolume(int)));
    ui->wordlineEdit->setText("Words: "+QString::number(size));
    qDebug()<<"size="<<this->size;

    timer=new QTimer();
    timer->start(1000); //timer will emit timeout() every second
    //connect(&timer, SIGNAL(timeout()), this, SLOT(slotTimeout()));
    QTime time_obj;
    //QTime time = ui->timeEdit->time().addSecs(-1);
    ui->timeEdit->setTime(time_obj.currentTime());
}

MainWindow::~MainWindow()
{
    delete ui;
    delete music;
}

//autodeFault on makes the buttonreact on the Enetre Key press
void MainWindow::on_pushButton_clicked()
{
    QString input=ui->lineEdit->text();
    qDebug()<<"input="<<input<<" QL size="<<names.size()<<" curnum="<<curnum;
    QString right_name=names[curnum];
    qDebug()<<"input="<<input<<" right_name="<<right_name
            <<" curnum="<<curnum<<" size="<<size;
    if(input == right_name){
        curnum=rand()%size;
        QString name=names[curnum];
        QPixmap pix=Map[name];
        ui->label->setPixmap(pix);
        ui->lineEdit->clear();
        qDebug()<<"A"<<" name="<<name;
    }
    else
    {
        QMessageBox* pmbx =
                            new QMessageBox("MessageBox",
                            "<b>Incorrect answer</b> <i>Simple</i>   <u>Message</u>",
                            QMessageBox::Warning,
                            QMessageBox::Yes,
                            QMessageBox::No,
                            QMessageBox::Cancel | QMessageBox::Escape);
        int n = pmbx->exec();
        delete pmbx;
        if (n == QMessageBox::Cancel || n == QMessageBox::Escape)
        {

        }
        qDebug()<<"B";
    }
}
//exit on pressing Escape
void MainWindow::keyPressEvent(QKeyEvent *e) {
    int key = e->key();
    if (key == Qt::Key_Escape) {
        exit(0);
    }
    QWidget::keyPressEvent(e);
}

void MainWindow::on_pushButton_2_clicked()
{
    //play backgound music
    //qDebug()<<" state="<<music->state();
    if(music->state() != QMediaPlayer::PlayingState)
    {
        music->play();
    }
}

void MainWindow::on_pushButton_3_clicked()
{
    music->pause();
//    if(music->state() == QMediaPlayer::PlayingState)
//    {
//        music->setPosition(0);
//    }
//    else if(music->state() == QMediaPlayer::StoppedState)
//        music->play();
}

void MainWindow::on_horizontalSlider_actionTriggered(int action)
{
    music->setVolume(action);
}

void MainWindow::on_quitButton_clicked()
{
    QMessageBox::StandardButton reply;
      reply = QMessageBox::question(this, "Quit?", "Quit?",
                                    QMessageBox::Yes|QMessageBox::No);
      if (reply == QMessageBox::Yes) {
        qDebug() << "Yes was clicked";
        QApplication::quit();
      } else {
        qDebug() << "Yes was *not* clicked";
      }
}


void MainWindow::on_actionNew_word_triggered()
{
    QMessageBox::information(this, "Required:", "The image should be 200x200, the format - .png or .jpg");
    //Создаем объект растрового изображения размером 200x200 пикселов (объект pix)
    QPixmap pix(200, 200);
    //Создаем объект строкового типа strFormat, в эту строку будет помещен
    //выбранный пользователем при помощи диалогового окна формат.
    QString strFilter="*.png";
    //имя файла:
    bool bOk;
    QString strname = QInputDialog::getText( 0,
                                         "Add a new word:",
                                         "Give the englich name to the picture:",
                                         QLineEdit::Normal,
                                         "name",
                                         &bOk
                                        );
    if (bOk) {
        //Последним передается адрес нашей строки, т. е. куда будет помещен
        //выбранный пользователем формат (объект strFilter)
        QString str = QFileDialog::getOpenFileName(0,
                                    "Сохранить изображение",
                                    "C:\\Users\\And\\Documents\\qtproj\\voc\\strname",
                                    "*.png ;; *.jpg ;; *.bmp",
                                    &strFilter
                                    );
        //После закрытия диалогового окна мы проверяем
        //строку str на содержимое, и если оно есть, то
        QPixmap pict;
        if (!str.isEmpty()) pict.load(str);
        Map.insert(strname,pict);
        names.push_back(strname);
        size=Map.size();
        ui->wordlineEdit->setText("Words: "+QString::number(size));
        qDebug()<<"str="<<str<<" strname="<<strname<<" size="<<size<<" bOk="<<bOk;
    }
    else QMessageBox::information(this, "info", "You did not give a name to the file");
}

void MainWindow::on_actionExit_triggered()
{
    QApplication::quit();
}
void MainWindow::slotTimeout()
{
    QTime time_obj;
    //QTime time = ui->timeEdit->time().addSecs(-1);
    ui->timeEdit->setTime(time_obj.currentTime());
    //if (time == QTime(0, 0))
        //time is zero, show message box
}

void MainWindow::on_timeEdit_userTimeChanged(const QTime &time)
{
    connect(timer, SIGNAL(timeout()), this, SLOT(slotTimeout()));
    timer->start(1000);
    slotTimeout();
}
