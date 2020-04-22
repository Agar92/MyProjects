#include "addressbookform.h"
#include "ui_addressbookform.h"

#include <QDesktopWidget>
#include <QSqlQueryModel>

AddressBookForm::AddressBookForm(QWidget *parent, DataBase _db) :
    QMainWindow(parent),
    ui(new Ui::AddressBookForm),
    db(_db)

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

AddressBookForm::~AddressBookForm()
{
    delete ui;
}

void AddressBookForm::SetText(QString name, QString surname)
{
    ui->labelPerson->setText(name + " " + surname);
}

void AddressBookForm::on_ExitButton_clicked()
{
    QApplication::quit();
}

void AddressBookForm::on_AddContactButton_clicked()
{
    //Dialog addContactDialog;
    addContactDialog = new Dialog();
    connect(addContactDialog, SIGNAL(DeleteDialog()), this, SLOT(CancelSlot()));
    connect(addContactDialog, SIGNAL(AddContact()), this, SLOT(AddContactSlot()));
    addContactDialog->setModal(true);
    addContactDialog->show();
    //addContactDialog.exec();
}

void AddressBookForm::CancelSlot()
{
    delete addContactDialog;
}

void AddressBookForm::AddContactSlot()
{
    QString name = addContactDialog->GetName();
    QString surname = addContactDialog->GetSurname();
    QString phone = addContactDialog->GetPhone();

    if(name.isEmpty() || surname.isEmpty() || phone.isEmpty() )
    {
        QMessageBox::critical(this, "Incomplete input", "Some of the input fields is empty");
    }
    else
    {
        addContactDialog->hide();
        delete addContactDialog;
        emit AddContactSignal(name, surname, phone);
    }
}

void AddressBookForm::on_ViewContactsButton_clicked()
{
    QSqlQueryModel * model = new QSqlQueryModel();
    QSqlQuery query(db.GetDB());
    query.exec("SELECT * FROM contacts;");
    QSqlRecord rec = query.record();
    int id;
    QString name, surname, phone;
    query.first();
    while (query.next()) {
        id = query.value(rec.indexOf("id")).toInt();
        name = query.value(rec.indexOf("name")).toString();
        surname = query.value(rec.indexOf("surname")).toString();
        phone = query.value(rec.indexOf("phone")).toString();


        qDebug() << "id is " << id
                 << ". name is " << name
                 << ". surname is " << surname
                 << " phone is " << phone
                 << query.lastError().text()
                 << query.lastError().number();
        }

    qDebug() << "STEP #3";

    model->setQuery(query);
    ui->tableView->setModel(model);

    ui->tableView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    ui->tableView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    ui->tableView->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    //renames the column name:
    //model->setHeaderData(1, Qt::Horizontal, "234");

}

void AddressBookForm::on_SearchButton_clicked()
{
    QSqlQueryModel * model = new QSqlQueryModel();
    QSqlQuery query(db.GetDB());

    QString name = ui->lineEditName->text();
    QString surname = ui->lineEditSurname->text();

    query.exec("SELECT * FROM contacts WHERE name LIKE ? AND surname LIKE ?;");
    query.bindValue(0, "%" + name + "%");
    query.bindValue(1, "%" + surname + "%");
    QSqlRecord rec = query.record();
    bool b = query.exec();
    int id;
    QString phone;
    query.first();
    while (query.next()) {
        id = query.value(rec.indexOf("id")).toInt();
        name = query.value(rec.indexOf("name")).toString();
        surname = query.value(rec.indexOf("surname")).toString();
        phone = query.value(rec.indexOf("phone")).toString();


        qDebug() << "id is " << id
                 << ". name is " << name
                 << ". surname is " << surname
                 << ". phone is " << phone
                 << query.lastError().text()
                 << query.lastError().number();
        }

    qDebug() << "STEP #3";

    model->setQuery(query);
    ui->tableView->setModel(model);

    ui->tableView->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    ui->tableView->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    ui->tableView->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    ui->tableView->verticalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    //renames the column name:
    //model->setHeaderData(1, Qt::Horizontal, "234");

}
