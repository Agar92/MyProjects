#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QMessageBox>
#include <QDesktopWidget>
#include <QPixmap>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
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

}

MainWindow::~MainWindow()
{
    delete ui;
    delete registrationForm;
    delete addressBookForm;
}

bool MainWindow::IsRegistered(QString login, QString password)
{
    return false;
}


void MainWindow::on_ExitButton_clicked()
{
    QApplication::quit();
}

void MainWindow::mousePressEvent(QMouseEvent *event) {
    m_nMouseClick_X_Coordinate = event->x();
    m_nMouseClick_Y_Coordinate = event->y();
}

void MainWindow::mouseMoveEvent(QMouseEvent *event) {
    move(event->globalX()-m_nMouseClick_X_Coordinate,event->globalY()-m_nMouseClick_Y_Coordinate);
}


void MainWindow::on_LoginButton_clicked()
{
    QString login = ui->lineEditLogin->text();
    QString password = ui->lineEditPassword->text();
    bool isRegistered = db.IsRegistered(login, password);
    if(isRegistered)
    {
        //QMessageBox::warning(this, "Correct log in", "Correct log in: You are already registered!");
        hide();
        QString name = db.GetName(login, password);
        QString surname = db.GetSurname(login, password);
        addressBookForm = new AddressBookForm(nullptr, db);
        connect(addressBookForm, SIGNAL(AddContactSignal(QString, QString, QString)), this,  SLOT(RegisterNewContact(QString, QString, QString)));
        addressBookForm->SetText(name, surname);
        addressBookForm->show();
    }
    else
    {
        QMessageBox::warning(this, "Not registered user tries to log in", "You need to register first");
    }
}

void MainWindow::on_RegisterButton_clicked()
{
    registrationForm = new RegistrationForm();
    hide();
    registrationForm->show();
    connect(registrationForm, SIGNAL(RegisterButtonPressed()), this, SLOT(Register()));
}

void MainWindow::Register()
{
    qDebug() << "BEGIN Register()";
    QString name = registrationForm->GetName();
    QString surname = registrationForm->GetSurname();
    QDate date = registrationForm->GetDate();
    QString occupation = registrationForm->GetOccupation();
    int spinValue = registrationForm->GetSpinBoxValue();
    QString login = registrationForm->GetLogin();
    QString password = registrationForm->GetPassword();

    qDebug() << "name=" << name
             << "surname=" << surname
             << "date="<<date.toString()
             << "occupation="<<occupation
             << " spinValue="<<spinValue;

    bool exist = db.FindNameAndSurnameAmongUsers(name, surname);
    if(exist)
    {
        QMessageBox::warning(this, "Already registered!", "You are already registered!");
    }
    else
    {
        db.RegisterCredentials(name, surname, date, login, password);

    }


    qDebug() << "END Register()";
}

void MainWindow::RegisterNewContact(QString name, QString surname, QString phone)
{
    qDebug() << "!!!RegisterNewContact() is called!!!";
    qDebug() << "name="<<name
             << "surname="<<surname
             << "phone="<<phone;

    db.AddNewContact(name, surname, phone);

}
