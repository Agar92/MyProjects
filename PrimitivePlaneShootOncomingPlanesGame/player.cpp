#include <QKeyEvent>
#include <QGraphicsScene>
#include "player.h"
#include "bullet.h"
#include "enemy.h"

#include <QDebug>
Player::Player(QGraphicsItem *parent) : QObject(), QGraphicsPixmapItem(parent)
{
    bulletsound = new QMediaPlayer();
    bulletsound->setMedia(QUrl("qrc:/sounds/shot.mp3"));
    setPixmap(QPixmap(":/images/plane1.png"));
}

void Player::keyPressEvent(QKeyEvent *event)
{
    //qDebug()<<"Player knows that you pressed the key";
    if(event->key() == Qt::Key_Left)
    {
        if(pos().x()>0)
        setPos(x()-10,y());
    }
    else if(event->key() == Qt::Key_Right)
    {
        if(pos().x()+100<800)
        setPos(x()+10,y());
    }
//    else if(event->key() == Qt::Key_Up)
//    {
//        setPos(x(),y()-10);
//    }
//    else if(event->key() == Qt::Key_Down)
//    {
//        setPos(x(),y()+10);
//    }
    else if(event->key() == Qt::Key_Space)
    {
        //create the bullet
        //qDebug()<<"MyRect knows that the bullet has been created";
        Bullet * bullet = new Bullet();
        bullet->setPos(x()+50,y());
        scene()->addItem(bullet);
        //play bullet sound
        if(bulletsound->state() == QMediaPlayer::PlayingState)
        {
            bulletsound->setPosition(0);
        }
        else if(bulletsound->state() == QMediaPlayer::StoppedState)
            bulletsound->play();
    }
}

void Player::spawn()
{
    //create an enemy
    Enemy * enemy = new Enemy();
    scene()->addItem(enemy);
}
