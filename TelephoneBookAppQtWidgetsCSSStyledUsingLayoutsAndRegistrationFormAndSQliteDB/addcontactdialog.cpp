#include "addcontactdialog.h"
#include "ui_addcontactdialog.h"

#include <QMessageBox>
#include <QDesktopWidget>
#include <QPixmap>
#include <QDebug>

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);

    //Скрываем элементы управления и рамку окна...
    this->setWindowFlags(Qt::FramelessWindowHint);
    //Включаем прозрачность главной формы...
    setAttribute(Qt::WA_TranslucentBackground );
    //center the frame:
    this->move(QApplication::desktop()->screen()->rect().center()-this->rect().center());

}

Dialog::~Dialog()
{
    delete ui;
}

QString Dialog::GetName(){return ui->lineEditName->text();}

QString Dialog::GetSurname(){return ui->lineEditSurname->text();}

QString Dialog::GetPhone(){return ui->lineEditPhone->text();}

void Dialog::on_CancelButton_clicked()
{
    hide();
    emit DeleteDialog();
}

void Dialog::mousePressEvent(QMouseEvent *event) {
    m_nMouseClick_X_Coordinate = event->x();
    m_nMouseClick_Y_Coordinate = event->y();
}

void Dialog::mouseMoveEvent(QMouseEvent *event) {
    move(event->globalX()-m_nMouseClick_X_Coordinate,event->globalY()-m_nMouseClick_Y_Coordinate);
}


void Dialog::on_AddContactButton_clicked()
{
    emit AddContact();
}
