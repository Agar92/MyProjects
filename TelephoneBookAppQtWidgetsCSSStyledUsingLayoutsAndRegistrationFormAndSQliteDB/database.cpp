#include "database.h"

DataBase::DataBase()
{
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE");
    db.setDatabaseName("mydb.sqlite");
    if (!db.open())
    {
        qDebug() << "***ERROR in DataBase(): Database not opened!";
    }
    query = QSqlQuery(db);
    // DDL query
    QString str = "CREATE TABLE IF NOT EXISTS credentials ("
                  "number INTEGER PRIMARY KEY NOT NULL, "
                  "name TEXT, "
                  "surname TEXT,"
                  "date_of_birth TEXT,"
                  "login TEXT,"
                  "password TEXT"
                  ");";
    bool b = query.exec(str);
    if (!b)
    {
        qDebug() << "Something went wrong! Table credentials not created! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "Table credentials is created b="<<b
                 << " query.isValid()="
                 << query.isValid();
    }

    CreateContactsTable();

}

bool DataBase::IsRegistered(QString login, QString password)
{
    qDebug() << "DataBase::IsRegistered() BEGIN:";
    query.clear();
    query.prepare("SELECT login, password FROM credentials WHERE login=? AND password=?;");
    query.bindValue(0, login);
    query.bindValue(1, password);
    bool b = query.exec();
    bool IsRegistered = false;
    if (!b)
    {
        qDebug() << "Something went wrong! SELECT credentials does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "SELECT credentials works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();

        query.first();

        qDebug() << query.value(0).toString()
                 << query.value(1).toString();

        if(query.value(0).toString() == login &&
                query.value(1).toString() == password)
            IsRegistered = true;

    }
    qDebug() << "DataBase::IsRegistered() END:";
    return IsRegistered;
}

bool DataBase::RegisterCredentials(QString name, QString surname,
                                   QDate birthdate, QString login,
                                   QString password)
{
    //query.clear();
    query.prepare("INSERT INTO credentials (name, surname, date_of_birth, login, password) "
                  "VALUES (?, ?, ?, ?, ?);");

    query.bindValue(0, name);
    query.bindValue(1, surname);
    query.bindValue(2, birthdate);
    query.bindValue(3, login);
    query.bindValue(4, password);
    bool b = query.exec();
    if (!b)
    {
        qDebug() << "Something went wrong! RegisterCredentials() does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
        return false;
    }
    else
    {
        qDebug() << "RegisterCredentials() works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();
        return true;

    }
}

bool DataBase::FindNameAndSurnameAmongUsers(QString name, QString surname)
{
    static int CALL=0;
    //check if already registered:
    query.clear();
    query.prepare("SELECT COUNT(*) FROM credentials WHERE name=? AND surname=?;");
    query.bindValue(0, name);
    query.bindValue(1, surname);
    bool b = query.exec();
    if (!b)
    {
        qDebug() << "Something went wrong! SELECT COUNT(*) in FindNameAndSurnameAmongUsers() does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
        return false;
    }
    else
    {
        qDebug() << "SELECT COUNT(*) in FindNameAndSurnameAmongUsers() works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();
    }
    int count=0;
    if (query.next()) {
        count = query.value(0).toInt();
    } else {
        qDebug() << query.lastError(); //check your error here
    }

    qDebug() << "count=" << count <<" CALL=" << CALL;
    ++CALL;

    if(count) return true;
    else      return false;
}

QString DataBase::GetName(QString login, QString password)
{
    query.clear();
    query.prepare("SELECT name FROM credentials WHERE login=? AND password=?;");
    query.bindValue(0, login);
    query.bindValue(1, password);
    bool b = query.exec();
    if (!b)
    {
        qDebug() << "Something went wrong! GetName() does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "GetName() works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();

    }
    query.first();
    qDebug() << query.value(0).toString();

    QString name = query.value(0).toString();
    return name;
}

QString DataBase::GetSurname(QString login, QString password)
{
    query.clear();
    query.prepare("SELECT surname FROM credentials WHERE login=? AND password=?;");
    query.bindValue(0, login);
    query.bindValue(1, password);
    bool b = query.exec();
    if (!b)
    {
        qDebug() << "Something went wrong! GetSurname() does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "GetSurname() works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();

    }
    query.first();
    qDebug() << query.value(0).toString();

    QString surname = query.value(0).toString();
    return surname;
}

void DataBase::CreateContactsTable()
{
    query.clear();
    QString str = "CREATE TABLE IF NOT EXISTS contacts ("
                  "id INTEGER PRIMARY KEY NOT NULL, "
                  "name TEXT, "
                  "surname TEXT,"
                  "phone TEXT"
                  ");";
    bool b = query.exec(str);
    if (!b)
    {
        qDebug() << "Something went wrong! Table contacts not created! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "Table contacts is created b="<<b
                 << " query.isValid()="
                 << query.isValid();
    }
}

void DataBase::AddNewContact(QString name, QString surname, QString phone)
{
    query.clear();
    query.prepare("INSERT INTO contacts (name, surname, phone) "
                  "VALUES (?, ?, ?);");

    query.bindValue(0, name);
    query.bindValue(1, surname);
    query.bindValue(2, phone);
    bool b = query.exec();
    if (!b)
    {
        qDebug() << "Something went wrong! AddNewContact() does not work! b="<<b
                 << query.lastError().isValid()
                 << query.lastError().text()
                 << query.lastError().number();
    }
    else
    {
        qDebug() << "AddNewContact() works properly b="<<b
                 << " query.isValid()="
                 << query.isValid();
    }
}
