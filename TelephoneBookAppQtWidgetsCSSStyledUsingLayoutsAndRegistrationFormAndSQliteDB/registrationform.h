#ifndef REGISTRATIONFORM_H
#define REGISTRATIONFORM_H

#include <QMainWindow>
#include <QMouseEvent>
#include <QDate>

namespace Ui {
class RegistrationForm;
}

class RegistrationForm : public QMainWindow
{
    Q_OBJECT

public:
    explicit RegistrationForm(QWidget *parent = nullptr);
    ~RegistrationForm();

    QString GetName();
    QString GetSurname();
    QDate GetDate();
    QString GetOccupation();
    int GetSpinBoxValue();
    QString GetLogin();
    QString GetPassword();


private slots:
    void on_ExitButton_clicked();

    void onregisterButtonclicked();

signals:
    void RegisterButtonPressed();

private:
    Ui::RegistrationForm *ui;

    //this is for moving the window with
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    int m_nMouseClick_X_Coordinate;
    int m_nMouseClick_Y_Coordinate;

};

#endif // REGISTRATIONFORM_H
