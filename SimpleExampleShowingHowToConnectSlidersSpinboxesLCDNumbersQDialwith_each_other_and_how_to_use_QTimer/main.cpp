#include <QtSql/QSql>
#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlRecord>
#include <QSqlDatabase>
#include <QSqlError>
#include <QDebug>
#include <QApplication>

#include "mainwindow.h"

int main(int argc, char *argv[])
{

    QApplication app(argc, argv);

    /*
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE");
        db.setDatabaseName("mydb.sqlite");
        if (!db.open()) {
            qDebug() << "Database not opened!";
            return -1;
        }

    QSqlQuery query(db);
    // DDL query
    QString str = "CREATE TABLE IF NOT EXISTS mytable ("
                "number INTEGER PRIMARY KEY NOT NULL, "
                "address TEXT, "
                "age INTEGER"
                ");";
    bool b = query.exec(str);
    if (!b) {
            qDebug() << "Something went wrong! b="<<b
                     << query.lastError().isValid()
                     << query.lastError().text()
                     << query.lastError().number();
    }
    else
        qDebug() << "Table is created b="<<b;

    QString str_insert = "INSERT INTO mytable(number, address, age)"
                "VALUES (%1, '%2', %3);";
    str = str_insert.arg("14").arg("hello world str.").arg("37");
    b = query.exec(str);
    if (!b) {
        qDebug() << "INSERT does not work!"
                 << query.lastError().text()
                 << query.lastError().number()
                 << "b="<<b;
    }
    if (!query.exec("SELECT * FROM mytable;")) {
        qDebug() << "Even SELECT does not work!"
                 << query.lastError().text()
                 << query.lastError().number()
                 << "b="<<b;
        return -2;
    }
    QSqlRecord rec = query.record();
    int number = 0, age = 0;
    QString address = "";

    while (query.next()) {
        number = query.value(rec.indexOf("number")).toInt();
        age = query.value(rec.indexOf("age")).toInt();
        address = query.value(rec.indexOf("address")).toString();

        qDebug() << "number is " << number
                 << ". age is " << age
                 << ". address" << address
                 << query.lastError().text()
                 << query.lastError().number();
        }

    db.close();
    */
    MainWindow w;
    w.show();

    return app.exec();
}
