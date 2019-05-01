#include "enemy.h"
#include "health.h"
#include "game.h"
#include <QTimer>
#include <QObject>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QDebug>
#include <stdlib.h>//rand() - large int.

extern Game * game;

Enemy::Enemy(QGraphicsItem * parent) : QObject(), QGraphicsPixmapItem(parent)
{
    //set random position
    int random_number=rand()%700;
    setPos(random_number,0);
    //draw the rect
    setPixmap(QPixmap(":/images/plane2.jpg"));
    setTransformOriginPoint(50,50);
    setRotation(180);
    //draw the rect
    //setRect(0,0,100,100);
    //connect
    QTimer * timer = new QTimer();
    connect(timer,SIGNAL(timeout()), this, SLOT(move()));
    timer->start(50);
}

void Enemy::move()
{
    //move enemy down
    setPos(x(),y()+5);
    if(pos().y()>600){
        game->health->decrease();
        scene()->removeItem(this);
        delete this;
        qDebug()<<"enemy deleted";
    }
}
