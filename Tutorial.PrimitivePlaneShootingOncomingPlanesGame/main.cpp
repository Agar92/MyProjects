#include <QApplication>
#include <QGraphicsScene>
#include "player.h"
#include <QGraphicsView>
#include <QTimer>
#include "game.h"

Game * game;

//the project is taken from: https://github.com/MeLikeyCode/QtGameTutorial/blob/master/tutorial8/
//YouTube channel: https://www.youtube.com/channel/UClzV7jGJREjvCTzfGTrdrkQ/videos
//Abdullah Aghazadah C++ Qt Game Tutorial 1 - Drawing the Player.
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    game = new Game();
    game->show();
    return a.exec();
}
