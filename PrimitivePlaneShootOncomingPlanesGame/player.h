#ifndef PLAYER_H
#define PLAYER_H
#include <QGraphicsPixmapItem>
#include <QObject>
#include <QMediaPlayer>
class Player : public QObject, public QGraphicsPixmapItem{
    Q_OBJECT
public:
    Player(QGraphicsItem * parent=0);
    void keyPressEvent(QKeyEvent * e);
public slots:
    void spawn();
private:
    QMediaPlayer * bulletsound;
};

#endif // PLAYER_H
