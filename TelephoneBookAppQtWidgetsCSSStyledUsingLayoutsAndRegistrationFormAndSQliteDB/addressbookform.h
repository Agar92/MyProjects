#ifndef ADDRESSBOOKFORM_H
#define ADDRESSBOOKFORM_H

#include <QMainWindow>

#include "addcontactdialog.h"
#include "database.h"

namespace Ui {
class AddressBookForm;
}

class AddressBookForm : public QMainWindow
{
    Q_OBJECT

public:
    explicit AddressBookForm(QWidget *parent = nullptr, DataBase _db = DataBase());
    ~AddressBookForm();
    void SetText(QString name, QString surname);

private slots:
    void on_ExitButton_clicked();

    void on_AddContactButton_clicked();

    void CancelSlot();

    void AddContactSlot();

    void on_ViewContactsButton_clicked();

    void on_SearchButton_clicked();

signals:
    void AddContactSignal(QString, QString, QString);

private:
    Ui::AddressBookForm *ui;

    Dialog * addContactDialog;

    QString name;
    QString surname;

    DataBase db;
};

#endif // ADDRESSBOOKFORM_H
