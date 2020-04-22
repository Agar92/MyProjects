#include "registrationform.h"
#include "ui_registrationform.h"

#include <QMessageBox>
#include <QDesktopWidget>
#include <QPixmap>

RegistrationForm::RegistrationForm(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::RegistrationForm)
{
    ui->setupUi(this);

    //Скрываем элементы управления и рамку окна...
    this->setWindowFlags(Qt::FramelessWindowHint);
    //Включаем прозрачность главной формы...
    setAttribute(Qt::WA_TranslucentBackground );
    //center the frame:
    this->move(QApplication::desktop()->screen()->rect().center()-this->rect().center());

    QPixmap pic(":/rec/img/close.png");
    int w = ui->ExitButton->width();
    int h = ui->ExitButton->height();
    pic.scaled(w, h, Qt::KeepAspectRatio);
    ui->ExitButton->setIcon(pic);

    //set QDateEdit format:
    QString format = QString("d MMMM yyyy");
    ui->dateEdit->setDisplayFormat(format);

    connect(ui->RegisterButton, SIGNAL(clicked()), this, SLOT(onregisterButtonclicked()));
}

RegistrationForm::~RegistrationForm()
{
    delete ui;
}

QString RegistrationForm::GetName(){return ui->lineEditName->text();}
QString RegistrationForm::GetSurname(){return ui->lineEditSurname->text();}
QDate RegistrationForm::GetDate(){return ui->dateEdit->date();}
QString RegistrationForm::GetOccupation(){return ui->comboBox->currentText();}
int RegistrationForm::GetSpinBoxValue(){return ui->spinBox->value();}
QString RegistrationForm::GetLogin(){return ui->lineEditLogin->text();}
QString RegistrationForm::GetPassword(){return ui->lineEditPassword->text();}


void RegistrationForm::on_ExitButton_clicked()
{
    QApplication::quit();
}

void RegistrationForm::mousePressEvent(QMouseEvent *event) {
    m_nMouseClick_X_Coordinate = event->x();
    m_nMouseClick_Y_Coordinate = event->y();
}

void RegistrationForm::mouseMoveEvent(QMouseEvent *event) {
    move(event->globalX()-m_nMouseClick_X_Coordinate,event->globalY()-m_nMouseClick_Y_Coordinate);
}

void RegistrationForm::onregisterButtonclicked()
{
    emit RegisterButtonPressed();
}
