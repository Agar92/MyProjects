#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <QMouseEvent>

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = nullptr);
    ~Dialog();

    QString GetName();
    QString GetSurname();
    QString GetPhone();

private slots:
    void on_CancelButton_clicked();

    void on_AddContactButton_clicked();

signals:
    void DeleteDialog();
    void AddContact();

private:
    Ui::Dialog *ui;

    //this is for moving the window with
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    int m_nMouseClick_X_Coordinate;
    int m_nMouseClick_Y_Coordinate;
};

#endif // DIALOG_H
