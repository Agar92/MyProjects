#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMouseEvent>

#include "registrationform.h"
#include "database.h"
#include "addressbookform.h"

#include <QtSql/QSqlTableModel>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    bool IsRegistered(QString login, QString password);

private slots:
    void on_ExitButton_clicked();

    void on_LoginButton_clicked();

    void on_RegisterButton_clicked();

    void Register();

    void RegisterNewContact(QString name, QString surname, QString phone);

private:
    Ui::MainWindow *ui;

    RegistrationForm * registrationForm;

    AddressBookForm * addressBookForm;

    DataBase db;

    //this is for moving the window with
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    int m_nMouseClick_X_Coordinate;
    int m_nMouseClick_Y_Coordinate;

};
#endif // MAINWINDOW_H
