#include "bullet.h"
#include "enemy.h"
#include <QTimer>
#include <QObject>
#include <QGraphicsScene>
#include <QDebug>
#include <QList>
#include "game.h"

extern Game * game;

Bullet::Bullet(QGraphicsItem *parent): QObject(), QGraphicsPixmapItem(parent)
{
    //draw the rect
    //setRect(0,0,10,50);
    setPixmap(QPixmap(":/images/bullet.png"));
    //connect
    QTimer * timer = new QTimer();
    connect(timer,SIGNAL(timeout()), this, SLOT(move()));
    timer->start(50);
}

void Bullet::move()
{
    //if bullet collides with the enemy, destroy both
    QList<QGraphicsItem *> colliding_items=collidingItems();
    for(int i=0,n=colliding_items.size(); i<n; ++i)
    {
        if(typeid(*(colliding_items[i])) == typeid(Enemy))
        {
            //increase the score
            game->score->increase();
            //remove them both
            scene()->removeItem(colliding_items[i]);
            scene()->removeItem(this);
            //delete thme both
            delete colliding_items[i];
            delete this;
            return;
        }
    }
    //move up
    setPos(x(),y()-10);
    //if(pos().y()+rect().height()<0){
    if(pos().y()<0){
        scene()->removeItem(this);
        delete this;
        qDebug()<<"bullet deleted";
    }
}
